"""Scale factor for 5P vs BODY subtraction (EndGraph).

Implements specs/scale-factor.md (approved 11/11/2026).
"""

import argparse
import json
import math
import sys
from bisect import bisect_right
from itertools import pairwise
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA_VERSION = 1
SKIP_PREFIXES = ("#", "track", "browser")


class ScaleFactorError(Exception):
    pass


@dataclass
class Segments:
    """Disjoint owned intervals of one (chrom, strand), sorted by start."""

    starts: list[int] = field(default_factory=list)
    ends: list[int] = field(default_factory=list)
    owners: list[tuple[str, ...]] = field(default_factory=list)


@dataclass
class Tally:
    """Sum terms for one track (5P or BODY)."""

    total: list[float] = field(default_factory=list)
    outside: list[float] = field(default_factory=list)
    # (gene, number of owners) -> terms; divided by the owner count once at the end
    shared: dict[tuple[str, int], list[float]] = field(default_factory=lambda: defaultdict(list))


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="degrseq scale-factor",
        description="Factor that 5P values are multiplied by before BODY is subtracted.",
    )
    parser.add_argument("--ends_plus", required=True, help="5P ends bedgraph, plus strand")
    parser.add_argument("--ends_minus", required=True, help="5P ends bedgraph, minus strand")
    parser.add_argument("--body_plus", required=True, help="BODY ends bedgraph, plus strand")
    parser.add_argument("--body_minus", required=True, help="BODY ends bedgraph, minus strand")
    parser.add_argument("--regions", required=True, help="gene regions, BED6")
    parser.add_argument("--fraglen", type=float, default=200)
    parser.add_argument("--body_strand", choices=("reverse", "forward"), default="reverse")
    parser.add_argument("--max_scale_factor", type=float, default=10)
    parser.add_argument("--write_genes", action="store_true", help="write per-gene table")
    parser.add_argument("--prefix", required=True, help="output file prefix")
    args = parser.parse_args(argv)
    if not (math.isfinite(args.fraglen) and args.fraglen > 0):
        parser.error("--fraglen must be > 0")
    if not (math.isfinite(args.max_scale_factor) and args.max_scale_factor > 0):
        parser.error("--max_scale_factor must be > 0")
    return args


def data_lines(path: str):
    """(line number, fields) for non-empty, non-header lines."""
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(SKIP_PREFIXES):
                continue
            yield number, line.split("\t")


def parse_interval(path: str, number: int, fields: list[str]) -> tuple[str, int, int]:
    try:
        start, end = int(fields[1]), int(fields[2])
    except (IndexError, ValueError):
        raise ScaleFactorError(f"{path}:{number}: bad coordinates") from None
    if start < 0 or end <= start:
        raise ScaleFactorError(f"{path}:{number}: start must be >= 0 and < end")
    return fields[0], start, end


def read_regions(path: str) -> tuple[dict[tuple[str, str], Segments], dict[str, int]]:
    """Owned segments per (chrom, strand) and gene lengths (rules 2, 5)."""
    rows: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    for number, fields in data_lines(path):
        if len(fields) < 6:
            raise ScaleFactorError(f"{path}:{number}: expected 6 columns (BED6)")
        chrom, start, end = parse_interval(path, number, fields)
        gene, strand = fields[3], fields[5]
        if strand not in ("+", "-"):
            raise ScaleFactorError(f"{path}:{number}: strand must be + or - (got {strand!r})")
        rows[(chrom, strand, gene)].append((start, end))

    lengths: dict[str, int] = defaultdict(int)
    events: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    for (chrom, strand, gene), intervals in rows.items():
        intervals.sort()
        for (_, prev_end), (start, _) in pairwise(intervals):
            if start < prev_end:
                raise ScaleFactorError(
                    f"{path}: rows of gene {gene} overlap on {chrom} {strand}; merge them first"
                )
        for start, end in intervals:
            lengths[gene] += end - start
            events[(chrom, strand)] += [(start, 1, gene), (end, -1, gene)]

    segments: dict[tuple[str, str], Segments] = {}
    for key, evs in events.items():
        evs.sort(key=lambda e: (e[0], e[1]))
        seg = Segments()
        active: set[str] = set()
        i = 0
        while i < len(evs):
            pos = evs[i][0]
            while i < len(evs) and evs[i][0] == pos:
                _, delta, gene = evs[i]
                (active.add if delta > 0 else active.discard)(gene)
                i += 1
            if active:
                seg.starts.append(pos)
                seg.ends.append(evs[i][0])
                seg.owners.append(tuple(sorted(active)))
        segments[key] = seg
    return segments, dict(lengths)


def tally_bedgraph(
    path: str,
    strand: str,
    segments: dict[tuple[str, str], Segments],
    tally: Tally,
    unknown_chroms: set[str],
    known_chroms: set[str],
) -> None:
    """Rules 2-4 for one bedgraph; a row adds its value at each nt it spans."""
    for number, fields in data_lines(path):
        if len(fields) < 4:
            raise ScaleFactorError(f"{path}:{number}: expected 4 columns")
        chrom, start, end = parse_interval(path, number, fields)
        try:
            value = float(fields[3])
        except ValueError:
            raise ScaleFactorError(f"{path}:{number}: bad value {fields[3]!r}") from None
        if not math.isfinite(value) or value < 0:
            raise ScaleFactorError(f"{path}:{number}: value must be finite and >= 0")
        if value == 0:
            continue
        if chrom not in known_chroms:
            unknown_chroms.add(chrom)
        owned = 0
        seg = segments.get((chrom, strand))
        if seg is not None:
            i = max(bisect_right(seg.starts, start) - 1, 0)
            while i < len(seg.starts) and seg.starts[i] < end:
                n = min(end, seg.ends[i]) - max(start, seg.starts[i])
                if n > 0:
                    owned += n
                    term = value * n
                    tally.total.append(term)
                    owners = seg.owners[i]
                    for gene in owners:
                        tally.shared[(gene, len(owners))].append(term)
                i += 1
        if owned < end - start:
            tally.outside.append(value * (end - start - owned))


def per_gene(tally: Tally) -> dict[str, float]:
    """Rule 4: shared terms summed per owner count, then divided once."""
    parts: dict[str, list[float]] = defaultdict(list)
    for (gene, k), terms in tally.shared.items():
        parts[gene].append(math.fsum(terms) / k)
    return {gene: math.fsum(p) for gene, p in parts.items()}


def compute(
    args: argparse.Namespace,
) -> tuple[dict, dict[str, int], dict[str, float], dict[str, float]]:
    segments, lengths = read_regions(args.regions)
    known = {chrom for chrom, _ in segments}
    unknown: set[str] = set()
    ends, body = Tally(), Tally()
    tally_bedgraph(args.ends_plus, "+", segments, ends, unknown, known)
    tally_bedgraph(args.ends_minus, "-", segments, ends, unknown, known)
    # rule 1: reverse = BODY aligned antisense to the gene
    body_for = {"+": args.body_minus, "-": args.body_plus}
    if args.body_strand == "forward":
        body_for = {"+": args.body_plus, "-": args.body_minus}
    for strand in ("+", "-"):
        tally_bedgraph(body_for[strand], strand, segments, body, unknown, known)

    end_total = math.fsum(ends.total)
    body_total = math.fsum(body.total)
    b = per_gene(body)
    with_body = {g: x for g, x in b.items() if x > 0}
    if end_total == 0:
        raise ScaleFactorError("no 5P values inside gene regions (E = 0)")
    if body_total == 0 or not with_body:
        raise ScaleFactorError("no BODY values inside gene regions (B = 0); check --body_strand")
    # rule 6
    mean_length = math.fsum(with_body.values()) / math.fsum(
        x / lengths[g] for g, x in with_body.items()
    )
    raw = (2 * args.fraglen / mean_length) * (body_total / end_total)
    factor = min(raw, args.max_scale_factor)
    factor_text = f"{factor:.4f}"
    stats = {
        "schema_version": SCHEMA_VERSION,
        "fraglen": args.fraglen,
        "body_strand": args.body_strand,
        "max_scale_factor": args.max_scale_factor,
        "end_total": end_total,
        "body_total": body_total,
        "mean_gene_length": mean_length,
        "genes_in_regions": len(lengths),
        "genes_with_body": len(with_body),
        "raw_factor": raw,
        "scale_factor": float(factor_text),
        "capped": raw > args.max_scale_factor,
        "values_outside_regions": {
            "end": math.fsum(ends.outside),
            "body": math.fsum(body.outside),
        },
        "chromosomes_without_regions": sorted(unknown),
    }
    return stats, lengths, per_gene(ends), b


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    try:
        stats, lengths, e, b = compute(args)
        Path(f"{args.prefix}.scale_factor.txt").write_text(
            f"{stats['scale_factor']:.4f}\n", encoding="utf-8"
        )
        if args.write_genes:
            with open(f"{args.prefix}.scale-factor.genes.tsv", "w", encoding="utf-8") as out:
                out.write("gene_id\tlength\tend\tbody\n")
                for gene in sorted(lengths, key=lambda g: g.encode()):
                    out.write(
                        f"{gene}\t{lengths[gene]}\t{e.get(gene, 0.0):.4f}\t{b.get(gene, 0.0):.4f}\n"
                    )
        Path(f"{args.prefix}.scale-factor.stats.json").write_text(
            json.dumps(stats, indent=2) + "\n", encoding="utf-8"
        )
    except (OSError, ScaleFactorError) as error:
        sys.exit(f"degrseq scale-factor: error: {error}")
