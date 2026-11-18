"""Capped / noncapped 5' end features of a replicate group (EndClass, numeric part).

Implements specs/cap-class.md (approved 18/11/2026; DL-074-DL-076, DL-080, DL-081, DL-083,
DL-085, DL-086).
"""

import argparse
import json
import math
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from degrseq.kde_peaks import KdePeaksError, data_lines, per_nt, read_bedgraph, read_lengths

SCHEMA_VERSION = 1
STRANDS = (("+", "plus"), ("-", "minus"))


class CapClassError(Exception):
    pass


@dataclass
class Region:
    chrom: str
    start: int
    end: int
    name: str
    strand: str


@dataclass
class Feature:
    chrom: str
    start: int  # trimmed
    end: int  # trimmed, exclusive
    strand: str
    name: str  # input region name
    peak: int  # offset from start
    total: float = 0.0
    ug: float = 0.0
    fraction: float = 0.0

    def key(self) -> tuple:
        return (self.chrom, self.start, self.end, 0 if self.strand == "+" else 1)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="degrseq cap-class",
        description="Sum replicates, trim features, split capped / noncapped by uG fraction.",
    )
    for flag, what in (
        ("--ends_plus", "5P ends (TSO-masked), plus"),
        ("--ends_minus", "5P ends (TSO-masked), minus"),
        ("--ug_plus", "cap (uG) track, plus"),
        ("--ug_minus", "cap (uG) track, minus"),
    ):
        parser.add_argument(
            flag, required=True, nargs="+", help=f"{what}: one bedgraph per replicate"
        )
    parser.add_argument("--features", required=True, help="replicable regions, BED6")
    parser.add_argument("--lengths", required=True, help="chromosome lengths (name, length)")
    parser.add_argument("--min_uug", type=float, default=0.1, help="capped if uG fraction >= this")
    parser.add_argument(
        "--min_reads", type=float, default=1.0, help="drop features with fewer 5P reads"
    )
    parser.add_argument("--decimals", type=int, default=4, help="decimals of the sum bedgraphs")
    parser.add_argument("--nanopare_compat", action="store_true", help="reproduce nanoPARE quirks")
    parser.add_argument("--prefix", required=True, help="output file prefix")
    args = parser.parse_args(argv)
    n = len(args.ends_plus)
    if not len(args.ends_minus) == len(args.ug_plus) == len(args.ug_minus) == n:
        parser.error(
            "--ends_plus, --ends_minus, --ug_plus, --ug_minus need one file per replicate each"
        )
    if not (math.isfinite(args.min_uug) and 0 <= args.min_uug <= 1):
        parser.error("--min_uug must be in [0, 1]")
    if not (math.isfinite(args.min_reads) and args.min_reads >= 0):
        parser.error("--min_reads must be >= 0")
    if not 0 <= args.decimals <= 10:
        parser.error("--decimals must be 0-10")
    return args


def read_regions(path: str, lengths: dict[str, int]) -> list[Region]:
    regions = []
    for number, fields in data_lines(path):
        try:
            chrom, start, end, name, strand = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                fields[3],
                fields[5],
            )
        except (IndexError, ValueError):
            raise CapClassError(f"{path}:{number}: expected BED6") from None
        if chrom not in lengths:
            raise CapClassError(f"{path}:{number}: chromosome {chrom!r} not in --lengths")
        if not 0 <= start < end <= lengths[chrom]:
            raise CapClassError(f"{path}:{number}: interval outside chromosome {chrom!r}")
        if strand not in ("+", "-"):
            raise CapClassError(f"{path}:{number}: strand must be + or -")
        regions.append(Region(chrom, start, end, name, strand))
    return regions


def round2_python(values: np.ndarray) -> np.ndarray:
    """Python round(x, 2) per value (exact decimal rounding)."""
    out = np.zeros_like(values)
    for i in np.flatnonzero(values):
        out[i] = round(float(values[i]), 2)
    return out


def read_tracks(paths: list[str], lengths: dict[str, int]) -> list[dict[str, list]]:
    return [read_bedgraph(p, lengths).rows for p in paths]


def chrom_sum(tracks: list[dict[str, list]], chrom: str, size: int, compat: bool) -> np.ndarray:
    """Rule 1: per-nt sum over replicates; compat: 2 decimals, position-0-only signal dropped."""
    arr = np.zeros(size)
    for rows in tracks:
        if chrom in rows:
            arr += per_nt(rows[chrom], size)
    if compat:
        arr = round2_python(arr)
        nz = np.flatnonzero(arr)
        if len(nz) == 1 and nz[0] == 0:  # quirk 37
            arr[0] = 0.0
    return arr


def bedgraph_rows(arr: np.ndarray, text) -> list[tuple[int, int, str]]:
    """Rule 2: non-zero positions, neighbours with equal written value merged."""
    rows: list[tuple[int, int, str]] = []
    for i in np.flatnonzero(arr):
        t = text(float(arr[i]))
        if rows and rows[-1][1] == i and rows[-1][2] == t:
            rows[-1] = (rows[-1][0], i + 1, t)
        else:
            rows.append((int(i), int(i) + 1, t))
    return rows


def auto_cutoff(fractions: list[float]) -> float:
    """Compat --min_uug 0: lowest-count 0.01 bin(s) among bins 0-24; mean of their lower edges."""
    hist = histogram(fractions)
    low = min(hist[:25])
    tied = [b for b in range(25) if hist[b] == low]
    return sum(b / 100 for b in tied) / len(tied)


def histogram(fractions: list[float]) -> list[int]:
    hist = [0] * 101
    for f in fractions:
        hist[min(100, math.floor(f * 100))] += 1
    return hist


def process(args: argparse.Namespace, lengths: dict[str, int]):
    compat = args.nanopare_compat
    text = (lambda v: str(v)) if compat else (lambda v: f"{v:.{args.decimals}f}")
    inputs = {
        ("5P", "+"): read_tracks(args.ends_plus, lengths),
        ("5P", "-"): read_tracks(args.ends_minus, lengths),
        ("uG", "+"): read_tracks(args.ug_plus, lengths),
        ("uG", "-"): read_tracks(args.ug_minus, lengths),
    }
    regions = read_regions(args.features, lengths)
    by_chrom: dict[str, list[Region]] = {}
    for r in regions:
        by_chrom.setdefault(r.chrom, []).append(r)
    stats = {"regions_in": len(regions), "dropped_zero_signal": 0, "dropped_below_min_reads": 0}
    sums: dict[tuple, dict[str, list]] = {k: {} for k in inputs}
    trimmed: list[Feature] = []
    kept: list[Feature] = []

    for chrom, size in lengths.items():  # one chromosome in memory at a time
        arr = {k: chrom_sum(v, chrom, size, compat) for k, v in inputs.items()}
        for k in inputs:
            sums[k][chrom] = bedgraph_rows(arr[k], text)
        feats = []
        for r in by_chrom.get(chrom, []):
            nz = np.flatnonzero(arr[("5P", r.strand)][r.start : r.end])
            if not len(nz):
                stats["dropped_zero_signal"] += 1
                continue
            a, b = r.start + int(nz[0]), r.start + int(nz[-1]) + 1
            span = arr[("5P", r.strand)][a:b]
            peak = int(len(span) - 1 - np.argmax(span[::-1]))  # rightmost maximum
            feats.append(Feature(chrom, a, b, r.strand, r.name, peak))
        feats.sort(key=Feature.key)
        trimmed += feats
        for f in feats:
            end = f.end + 1 if compat else f.end  # quirk 15
            five = arr[("5P", f.strand)][f.start : end]
            ug = arr[("uG", f.strand)][f.start : end]
            if compat:
                f.total, f.ug = sum(float(x) for x in five), sum(float(x) for x in ug)
            else:
                f.total, f.ug = math.fsum(five), math.fsum(ug)
            if f.total < args.min_reads:
                stats["dropped_below_min_reads"] += 1
                continue
            f.fraction = f.ug / f.total if f.total else 0.0
            if not compat and f.fraction > 1 + 1e-9:
                raise CapClassError(
                    f"uG exceeds 5P in feature {chrom}:{f.start}-{f.end} {f.strand}; "
                    "cap track counted per G (compat end-counts)? "
                    "use --nanopare_compat for both (DL-074)"
                )
            kept.append(f)

    cutoff = args.min_uug
    if compat and cutoff == 0:
        cutoff = auto_cutoff([f.fraction for f in kept])
    capped = [f for f in kept if f.fraction >= cutoff]
    noncapped = [f for f in kept if f.fraction < cutoff]
    stats.update(
        schema_version=SCHEMA_VERSION,
        mode="nanopare_compat" if compat else "default",
        replicates=len(args.ends_plus),
        min_uug=cutoff,
        min_reads=args.min_reads,
        capped={s: sum(f.strand == s for f in capped) for s, _ in STRANDS},
        noncapped={s: sum(f.strand == s for f in noncapped) for s, _ in STRANDS},
        ug_fraction_histogram=histogram([f.fraction for f in kept]),
    )
    return sums, trimmed, capped, noncapped, stats


def write_outputs(args, lengths, sums, trimmed, capped, noncapped, stats) -> None:
    compat = args.nanopare_compat
    chroms = sorted(lengths) if compat else list(lengths)
    for kind, tag in (("5P", ""), ("uG", "uG.")):
        for st, label in STRANDS:
            with open(f"{args.prefix}.{tag}{label}.bedgraph", "w", encoding="utf-8") as out:
                for chrom in chroms:
                    for a, b, t in sums[(kind, st)][chrom]:
                        out.write(f"{chrom}\t{a}\t{b}\t{t}\n")
    with open(f"{args.prefix}.peaks.bed", "w", encoding="utf-8") as out:
        for f in trimmed:
            out.write(f"{f.chrom}\t{f.start}\t{f.end}\t{f.name}\t{f.peak}\t{f.strand}\n")
    for label, feats in (("capped", capped), ("noncapped", noncapped)):
        with open(f"{args.prefix}.{label}.bed", "w", encoding="utf-8") as out:
            for n, f in enumerate(feats, 1):
                end, peak = (f.end + 1, f.peak + 1) if compat else (f.end, f.peak)
                frac = str(round(f.fraction, 3)) if compat else f"{f.fraction:.3f}"
                out.write(
                    f"{f.chrom}\t{f.start}\t{end}\t5P.{label}.{n}\t{peak}\t{f.strand}\t.\t.\t.\t{frac}\n"
                )
    Path(f"{args.prefix}.cap-class.stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    try:
        lengths = read_lengths(args.lengths)
        sums, trimmed, capped, noncapped, stats = process(args, lengths)
    except (CapClassError, KdePeaksError, OSError) as error:
        print(f"degrseq cap-class: error: {error}", file=sys.stderr)
        sys.exit(1)
    write_outputs(args, lengths, sums, trimmed, capped, noncapped, stats)
