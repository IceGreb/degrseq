"""5' end features by subtractive kernel density estimation (EndGraph).

Implements specs/kde-peaks.md (approved 16/11/2026, DL-072).
"""

import argparse
import json
import math
import sys
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

SCHEMA_VERSION = 1
SKIP_PREFIXES = ("#", "track", "browser")
STRANDS = (("+", "plus"), ("-", "minus"))


class KdePeaksError(Exception):
    pass


@dataclass
class Feature:
    chrom: str
    start: int
    end: int  # exclusive
    strand: str
    name: str
    dominant: str
    rpm: float
    rpm_text: str
    peak: int
    other_peaks: str


@dataclass
class StrandStats:
    positive_positions: int = 0
    features_called: int = 0
    features_kept: int = 0
    features_without_reads: int = 0


@dataclass
class Bedgraph:
    """Rows per chromosome, in file order."""

    rows: dict[str, list[tuple[int, int, float]]] = field(default_factory=lambda: defaultdict(list))


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="degrseq kde-peaks",
        description="5' end features from scale x 5P - BODY by Laplace kernel density.",
    )
    parser.add_argument("--ends_plus", required=True, help="5P ends bedgraph (TSO-masked), plus")
    parser.add_argument("--ends_minus", required=True, help="5P ends bedgraph (TSO-masked), minus")
    parser.add_argument("--body_plus", required=True, help="BODY ends bedgraph, plus")
    parser.add_argument("--body_minus", required=True, help="BODY ends bedgraph, minus")
    parser.add_argument("--scale_factor", required=True, help="file with the scale factor")
    parser.add_argument("--lengths", required=True, help="chromosome lengths (name, length)")
    parser.add_argument("--bandwidth", type=int, default=15)
    parser.add_argument("--kernel_extent", type=int, default=3)
    parser.add_argument("--decimals", type=int, default=3)
    parser.add_argument("--min_length", type=int, default=10)
    parser.add_argument("--min_rpm", type=float, default=0.5)
    parser.add_argument(
        "--body_strand", choices=("reverse", "forward", "unstranded"), default="reverse"
    )
    parser.add_argument("--write_density", action="store_true", help="write density bedgraphs")
    parser.add_argument("--nanopare_compat", action="store_true", help="reproduce nanoPARE quirks")
    parser.add_argument("--prefix", required=True, help="output file prefix")
    args = parser.parse_args(argv)
    if not 1 <= args.bandwidth <= 30:
        parser.error("--bandwidth must be 1-30")
    if not 1 <= args.kernel_extent <= 6:
        parser.error("--kernel_extent must be 1-6")
    if not 0 <= args.decimals <= 10:
        parser.error("--decimals must be 0-10")
    if args.min_length < 1:
        parser.error("--min_length must be >= 1")
    if not (math.isfinite(args.min_rpm) and args.min_rpm >= 0):
        parser.error("--min_rpm must be >= 0")
    return args


def data_lines(path: str):
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(SKIP_PREFIXES):
                continue
            yield number, line.split("\t")


def read_lengths(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for number, fields in data_lines(path):
        try:
            name, length = fields[0], int(fields[1])
        except (IndexError, ValueError):
            raise KdePeaksError(f"{path}:{number}: expected name<TAB>length") from None
        if length < 1 or name in lengths:
            raise KdePeaksError(f"{path}:{number}: invalid or duplicate chromosome {name!r}")
        lengths[name] = length
    return lengths


def read_bedgraph(path: str, lengths: dict[str, int]) -> Bedgraph:
    bg = Bedgraph()
    for number, fields in data_lines(path):
        try:
            chrom, start, end, value = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
        except (IndexError, ValueError):
            raise KdePeaksError(f"{path}:{number}: malformed bedgraph row") from None
        if chrom not in lengths:
            raise KdePeaksError(f"{path}:{number}: chromosome {chrom!r} not in --lengths")
        if not 0 <= start < end <= lengths[chrom]:
            raise KdePeaksError(f"{path}:{number}: interval outside chromosome {chrom!r}")
        if not (math.isfinite(value) and value >= 0):
            raise KdePeaksError(f"{path}:{number}: value must be finite and >= 0")
        if value:
            bg.rows[chrom].append((start, end, value))
    return bg


def read_scale_factor(path: str) -> float:
    try:
        value = float(Path(path).read_text(encoding="utf-8").split()[0])
    except (IndexError, ValueError, OSError):
        raise KdePeaksError(f"{path}: cannot read a scale factor") from None
    if not (math.isfinite(value) and value > 0):
        raise KdePeaksError(f"{path}: scale factor must be finite and > 0")
    return value


def per_nt(rows: list[tuple[int, int, float]], size: int) -> np.ndarray:
    arr = np.zeros(size)
    for start, end, value in rows:
        arr[start:end] += value
    return arr


def round_half_even(values: np.ndarray, decimals: int) -> np.ndarray:
    """Round as Python's round() (exact decimal); np.round only differs next to a half."""
    out = np.round(values, decimals)
    scaled = values * 10.0**decimals
    near = np.flatnonzero(np.abs(scaled - np.floor(scaled) - 0.5) < 1e-6)
    for i in near:
        out[i] = round(float(values[i]), decimals)
    return out


def kernel(bandwidth: int, extent: int, compat: bool) -> np.ndarray:
    """Rule 3: Laplace weights for offsets -W..W, normalised; compat: single precision, drops +W."""
    w = bandwidth * extent
    if compat:
        # DL-071: exp evaluated in float32, then /2H and normalised in float64
        x = np.arange(-w, w + 1, dtype=np.float32)
        k = np.exp(-(np.abs(x) / np.float32(bandwidth))).astype(np.float64) / (2 * bandwidth)
        return (k / k.sum())[:-1]
    x = np.arange(-w, w + 1, dtype=np.float64)
    k = np.exp(-np.abs(x) / bandwidth)
    return k / k.sum()


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """[start, end) of consecutive True values."""
    edges = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    return list(
        zip(np.flatnonzero(edges == 1).tolist(), np.flatnonzero(edges == -1).tolist(), strict=True)
    )


def density_peaks(values: Sequence[float]) -> list[int]:
    """Rule 7: last position of each equal-value run entered by a rise and followed by a fall."""
    groups: list[list] = []
    for i, v in enumerate(values):
        if groups and groups[-1][2] == v:
            groups[-1][1] = i
        else:
            groups.append([i, i, v])
    peaks, previous = [], 0.0
    rises = []
    for _, _, v in groups:
        rises.append(v > previous)
        previous = v
    for j, (_, end, v) in enumerate(groups):
        if rises[j] and j + 1 < len(groups) and groups[j + 1][2] < v:
            peaks.append(end)
    if not peaks and groups and rises[-1]:
        peaks = [len(values) - 1]
    return sorted(peaks, key=lambda p: (values[p], p), reverse=True)


def build_difference(
    ends: np.ndarray, body: np.ndarray, scale: float, compat: bool, covered: np.ndarray
) -> np.ndarray:
    """Rule 2: S x E - B; compat rounds per nt to 2 decimals and drops a position-0-only track."""
    d = ends * scale - body
    if compat:
        nz = np.flatnonzero(d)
        d[nz] = round_half_even(d[nz], 2)
        cov = np.flatnonzero(covered)
        if cov.size == 1 and cov[0] == 0:
            d[0] = 0.0
    return d


def process(args: argparse.Namespace) -> tuple[list[Feature], dict, dict]:
    compat = args.nanopare_compat
    lengths = read_lengths(args.lengths)
    scale = read_scale_factor(args.scale_factor)
    ends = {
        "+": read_bedgraph(args.ends_plus, lengths),
        "-": read_bedgraph(args.ends_minus, lengths),
    }
    body = {
        "+": read_bedgraph(args.body_plus, lengths),
        "-": read_bedgraph(args.body_minus, lengths),
    }
    # rule 1
    pair = {"+": ["-"], "-": ["+"]}
    if args.body_strand == "forward":
        pair = {"+": ["+"], "-": ["-"]}
    elif args.body_strand == "unstranded":
        pair = {"+": ["+", "-"], "-": ["+", "-"]}

    # rule 9: total 5P
    if compat:
        plus_total = 0.0
        for chrom in lengths:
            for start, end, value in ends["+"].rows.get(chrom, []):
                for _ in range(end - start):
                    plus_total += value
        total = float(int(plus_total + plus_total))
    else:
        total = math.fsum(
            v * (e - s) for st in "+-" for rows in ends[st].rows.values() for s, e, v in rows
        )

    k = kernel(args.bandwidth, args.kernel_extent, compat)
    w = args.bandwidth * args.kernel_extent
    features: list[Feature] = []
    density: dict[str, list[tuple[str, int, int, float]]] = {"+": [], "-": []}
    stats = {st: StrandStats() for st, _ in STRANDS}

    for chrom, length in lengths.items():
        for st, label in STRANDS:
            e_rows = ends[st].rows.get(chrom, [])
            b_rows = [r for p in pair[st] for r in body[p].rows.get(chrom, [])]
            if not e_rows:
                continue  # no positive density without 5P ends
            size = length + w if compat else length
            e_arr = per_nt(e_rows, size)
            b_arr = per_nt(b_rows, size)
            covered = (e_arr != 0) | (b_arr != 0) if compat else None
            d = build_difference(e_arr[:length], b_arr[:length], scale, compat, covered)
            full = np.convolve(d, k)  # rule 4
            dens = full[w : w + size]
            dens = round_half_even(dens, args.decimals)  # rule 5
            positive = dens > 0
            stats[st].positive_positions += int(positive.sum())
            n = 0
            for start, end in runs(positive):
                if end - start < args.min_length:
                    continue
                n += 1
                stats[st].features_called += 1
                vals = dens[start:end].tolist()
                peaks = density_peaks(vals)
                last = end - 1 if compat else end  # rule 8; compat leaves out the last nt
                window = e_arr[start:last]
                dom = str(int(window.size - 1 - np.argmax(window[::-1]))) if window.size else "."
                if compat:
                    reads = sum(window.tolist())
                    if total == 0:
                        raise KdePeaksError("total of --ends_plus is 0 (nanoPARE divides by zero)")
                    rpm = round(reads / total * 1e6, 2)
                    rpm_text = str(rpm)
                else:
                    reads = math.fsum(window.tolist())
                    rpm = reads / total * 1e6 if total else 0.0
                    rpm_text = f"{rpm:.2f}"
                if reads == 0:
                    stats[st].features_without_reads += 1
                if rpm < args.min_rpm:  # rule 10
                    continue
                stats[st].features_kept += 1
                features.append(
                    Feature(
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=st,
                        name=f"5P.{label}.{chrom}.{n}",
                        dominant=dom,
                        rpm=rpm,
                        rpm_text=rpm_text,
                        peak=peaks[0],
                        other_peaks=",".join(str(p) for p in peaks[1:10]) or "-",
                    )
                )
            if args.write_density:
                for start, end in runs(positive):
                    j = start
                    while j < end:
                        v = dens[j]
                        m = j + 1
                        while m < end and dens[m] == v:
                            m += 1
                        density[st].append((chrom, j, m, float(v)))
                        j = m

    order = {c: i for i, c in enumerate(lengths)}
    features.sort(key=lambda f: (order[f.chrom], f.start, f.end, f.strand != "+"))  # rule 11
    if compat:
        for st in density:
            density[st].sort(key=lambda r: (r[0], r[1]))  # original: chromosome names sorted
    out_stats = {
        "schema_version": SCHEMA_VERSION,
        "mode": "nanopare_compat" if compat else "default",
        "scale_factor": scale,
        "body_strand": args.body_strand,
        "bandwidth": args.bandwidth,
        "kernel_extent": args.kernel_extent,
        "decimals": args.decimals,
        "min_length": args.min_length,
        "min_rpm": args.min_rpm,
        "total_5p": total,
        "plus": vars(stats["+"]),
        "minus": vars(stats["-"]),
    }
    return features, density, out_stats


def write_outputs(
    args: argparse.Namespace, features: list[Feature], density: dict, stats: dict
) -> None:
    compat = args.nanopare_compat
    with open(f"{args.prefix}.end_features.bed", "w", encoding="utf-8") as out:
        for f in features:
            end = f.end - 1 if compat else f.end  # quirk 12
            out.write(
                f"{f.chrom}\t{f.start}\t{end}\t{f.name}\t{f.dominant}\t{f.strand}\t"
                f"{f.rpm_text}\t{f.peak}\t{f.other_peaks}\n"
            )
    if args.write_density:
        for st, label in STRANDS:
            with open(f"{args.prefix}.density_{label}.bedgraph", "w", encoding="utf-8") as out:
                for chrom, start, end, v in density[st]:
                    text = str(v) if compat else f"{v:.{args.decimals}f}"
                    out.write(f"{chrom}\t{start}\t{end}\t{text}\n")
    Path(f"{args.prefix}.kde-peaks.stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    try:
        features, density, stats = process(args)
    except KdePeaksError as error:
        print(f"degrseq kde-peaks: error: {error}", file=sys.stderr)
        sys.exit(1)
    write_outputs(args, features, density, stats)
