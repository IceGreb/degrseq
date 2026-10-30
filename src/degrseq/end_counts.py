"""Per-nucleotide RNA 5' end counts from STAR alignments.

Implements specs/end-counts.md (approved 18/10/2026): default mode and `--nanopare_compat`.
"""

import argparse
import json
import sys
from array import array
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam

SCHEMA_VERSION = 1

# Rule 2
SKIP_FLAGS = 0x4 | 0x200 | 0x400 | 0x800

# Rule 5: M, =, X consume reference and add positions; D, N consume reference only.
ALIGNED_OPS = {0, 7, 8}
REF_ONLY_OPS = {2, 3}
SOFT_CLIP, HARD_CLIP = 4, 5

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

READ_FATES = (
    "unique",
    "unique_after_filter",
    "multimapped_rescued",
    "multimapped_naive",
    "multimapped_uniform",
    "multimapped_dropped",
    "multimapped_lost_compat",
    "discarded_too_short",
    "discarded_untemplated_non_g",
    "discarded_3prime_clip_only",
    "discarded_no_mate1",
    "discarded_flagged",
    "unmapped",
)
LOCUS_COUNTERS = (
    "loci_too_short",
    "loci_untemplated_non_g",
    "loci_no_mate1",
    "loci_key_collisions",
)
G_RUN_BINS = ("0", "1", "2", "3", ">=4")
COVERAGE_CHUNK = 1_000_000  # intervals per vectorised add (compat coverage)


class EndCountsError(Exception):
    """Input that violates the spec's preconditions."""


@dataclass
class Locus:
    ref_id: int
    strand: str
    five: int
    intervals: list[tuple[int, int]]
    n_positions: int
    aln_start: int
    untemplated5: str = ""
    has_untemplated3: bool = False
    status: str = "ok"
    # Compat only: CIGAR positions, untemplated (coordinate, read-sense base) per side.
    n_cigar: int = 0
    cigar_intervals: tuple = ()
    five_sites: tuple = ()
    three_sites: tuple = ()
    cigar: str = ""
    secondary: bool = False

    @property
    def cap(self) -> bool:
        return len(self.untemplated5) > 0


@dataclass
class Block:
    loci: list[Locus]
    sort_key: tuple
    weight: float = 0.0
    n_reads: int = 0
    shares_idx: list[np.ndarray] = field(default_factory=list)  # locus positions in the grid


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="degrseq end-counts",
        description="Count RNA 5' ends per nucleotide from STAR alignments.",
    )
    parser.add_argument("--bam", required=True, help="name-grouped BAM (STAR Unsorted)")
    parser.add_argument("--fasta", help="genome FASTA with .fai (required for 5P)")
    parser.add_argument("--library", required=True, choices=("5P", "BODY"))
    parser.add_argument("--prefix", required=True, help="output file prefix")
    parser.add_argument("--min_aligned", type=int, default=16)
    parser.add_argument(
        "--multimapper_mode", choices=("rescue", "uniform", "unique"), default="rescue"
    )
    parser.add_argument(
        "--decimals", type=int, help="default 4; not allowed with --nanopare_compat"
    )
    parser.add_argument("--write_coverage", action="store_true")
    parser.add_argument("--nanopare_compat", action="store_true")
    args = parser.parse_args(argv)
    if args.nanopare_compat and args.multimapper_mode != "rescue":
        parser.error("--nanopare_compat requires --multimapper_mode rescue")
    if args.nanopare_compat and args.decimals is not None:
        parser.error("--decimals can't be used with --nanopare_compat (fixed at 2)")
    if args.decimals is None:
        args.decimals = 2 if args.nanopare_compat else 4
    if args.library == "5P" and args.fasta is None:
        parser.error("--fasta is required for --library 5P")
    if args.min_aligned < 1 or args.decimals < 0:
        parser.error("--min_aligned must be >= 1 and --decimals >= 0")
    return args


def read_groups(bam: pysam.AlignmentFile) -> Iterator[list[pysam.AlignedSegment]]:
    """Rule 1: consecutive records with the same query name."""
    group: list[pysam.AlignedSegment] = []
    for record in bam.fetch(until_eof=True):
        if group and record.query_name != group[0].query_name:
            yield group
            group = []
        group.append(record)
    if group:
        yield group


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


class Genome:
    def __init__(self, fasta: pysam.FastaFile | None, lengths: dict[str, int]) -> None:
        self.fasta = fasta
        self.lengths = lengths

    def fetch(self, ref: str, start: int, end: int) -> str:
        assert self.fasta is not None
        return self.fasta.fetch(ref, start, end).upper()


def record_locus(record: pysam.AlignedSegment, genome: Genome, reclaim: bool) -> Locus:
    """Rules 5-7 for one record: positions, reclaim, 5' end and untemplated bases."""
    ops = [(op, n) for op, n in record.cigartuples if op != HARD_CLIP]
    pos = record.reference_start
    intervals: list[tuple[int, int]] = []
    for op, n in ops:
        if op in ALIGNED_OPS:
            if intervals and intervals[-1][1] == pos:
                intervals[-1] = (intervals[-1][0], pos + n)
            else:
                intervals.append((pos, pos + n))
            pos += n
        elif op in REF_ONLY_OPS:
            pos += n
    if not intervals:
        raise EndCountsError(f"read {record.query_name}: mapped record with no aligned bases")
    start, end = intervals[0][0], intervals[-1][1]
    cigar_intervals = tuple(intervals)
    n_cigar = sum(e - s for s, e in intervals)

    left = ops[0][1] if ops[0][0] == SOFT_CLIP else 0
    right = ops[-1][1] if len(ops) > 1 and ops[-1][0] == SOFT_CLIP else 0
    left_untemp = right_untemp = ""
    if reclaim and (left or right):
        seq = record.query_sequence
        if seq is None:
            raise EndCountsError(f"read {record.query_name}: soft-clipped record without SEQ")
        seq = seq.upper()
        ref = record.reference_name
        if left:
            clip = seq[:left]
            genome_left = genome.fetch(ref, max(0, start - left), start)
            k = 0
            # Walk outward from the alignment edge; N never matches.
            while k < len(genome_left) and clip[-1 - k] == genome_left[-1 - k] != "N":
                k += 1
            left_untemp = clip[: left - k]
            start -= k
            if k:
                intervals[0] = (start, intervals[0][1])
        if right:
            clip = seq[len(seq) - right :]
            genome_right = genome.fetch(ref, end, min(genome.lengths[ref], end + right))
            k = 0
            while k < len(genome_right) and clip[k] == genome_right[k] != "N":
                k += 1
            right_untemp = clip[k:]
            end += k
            if k:
                intervals[-1] = (intervals[-1][0], end)

    left_sites = tuple(zip(range(start - len(left_untemp), start), left_untemp, strict=True))
    right_sites = tuple(zip(range(end, end + len(right_untemp)), right_untemp, strict=True))
    if record.is_reverse:
        strand, five = "-", end - 1
        untemp5 = right_untemp.translate(COMPLEMENT)[::-1]
        untemp3 = left_untemp
        five_sites = tuple((p, b.translate(COMPLEMENT)) for p, b in right_sites)
        three_sites = tuple((p, b.translate(COMPLEMENT)) for p, b in left_sites)
    else:
        strand, five = "+", start
        untemp5, untemp3 = left_untemp, right_untemp
        five_sites, three_sites = left_sites, right_sites
    return Locus(
        ref_id=record.reference_id,
        strand=strand,
        five=five,
        intervals=intervals,
        n_positions=sum(e - s for s, e in intervals),
        aln_start=record.reference_start,
        untemplated5=untemp5,
        has_untemplated3=bool(untemp3),
        n_cigar=n_cigar,
        cigar_intervals=cigar_intervals,
        five_sites=five_sites,
        three_sites=three_sites,
        cigar=record.cigarstring or "",
        secondary=record.is_secondary,
    )


def pair_records(
    mate1: list[pysam.AlignedSegment], mate2: list[pysam.AlignedSegment]
) -> list[tuple[pysam.AlignedSegment | None, pysam.AlignedSegment | None]]:
    """Rule 17: a locus is a mate-1 record plus the mate-2 record with the same HI tag."""

    def hit(record: pysam.AlignedSegment) -> tuple:
        if record.has_tag("HI"):
            return ("HI", record.get_tag("HI"))
        return ("pos", record.reference_id, record.reference_start)

    def mate_hit(record: pysam.AlignedSegment) -> tuple:
        if record.has_tag("HI"):
            return ("HI", record.get_tag("HI"))
        return ("pos", record.next_reference_id, record.next_reference_start)

    by_hit = {}
    for record in mate2:
        by_hit.setdefault(hit(record), record)
    pairs: list[tuple[pysam.AlignedSegment | None, pysam.AlignedSegment | None]] = []
    used = set()
    for record in mate1:
        partner = by_hit.get(mate_hit(record))
        if partner is not None:
            used.add(id(partner))
        pairs.append((record, partner))
    pairs.extend((None, record) for record in mate2 if id(record) not in used)
    return pairs


def pair_locus(
    m1: pysam.AlignedSegment | None, m2: pysam.AlignedSegment | None, genome: Genome, reclaim: bool
) -> Locus:
    if m1 is None:
        assert m2 is not None
        locus = record_locus(m2, genome, reclaim)
        locus.status = "no_mate1"
        return locus
    locus = record_locus(m1, genome, reclaim)
    if m2 is not None and m2.reference_id == m1.reference_id:
        # Mate 2 adds coverage on mate 1's strand; its own clips decide nothing.
        other = record_locus(m2, genome, reclaim)
        locus.intervals = merge_intervals(locus.intervals + other.intervals)
        locus.n_positions = sum(e - s for s, e in locus.intervals)
        locus.aln_start = min(locus.aln_start, other.aln_start)
    return locus


def read_weight(records: list[pysam.AlignedSegment]) -> float:
    """Rule 13."""
    values = {float(r.get_tag("XW")) for r in records if r.has_tag("XW")}
    if len(values) > 1:
        raise EndCountsError(f"read {records[0].query_name}: records with different XW values")
    return values.pop() if values else 1.0


def read_sequence(records: list[pysam.AlignedSegment]) -> str:
    """Mate 1 + mate 2 in read orientation, for block identity (rule 20)."""
    parts = []
    for want_read2 in (False, True):
        for r in records:
            if r.is_read2 == want_read2 and r.query_sequence is not None:
                parts.append(r.get_forward_sequence())
                break
    return "".join(parts)


def g_run_bin(n: int) -> str:
    return ">=4" if n >= 4 else str(n)


def format_value(value: float, decimals: int) -> str:
    text = f"{value:.{decimals}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


class Counter:
    def __init__(self, args: argparse.Namespace, refs: list[str], genome: Genome) -> None:
        self.args = args
        self.refs = refs
        self.genome = genome
        self.is_5p = args.library == "5P"
        self.compat = args.nanopare_compat
        self.ends: dict[tuple[int, str, int], float] = {}
        self.cap: dict[tuple[int, str, int], float] = {}
        self.stats: dict[str, float | int] = dict.fromkeys(
            ("reads_total", *READ_FATES, *LOCUS_COUNTERS, "kept_3prime_clip_only", "nh_mismatch"),
            0,
        )
        self.cap_positive_reads = 0.0
        self.g_runs = dict.fromkeys(G_RUN_BINS, 0)
        # Unique-read coverage, kept as intervals: ref, strand (0 +, 1 -), start, end, weight.
        self.cov_ref, self.cov_strand = array("i"), array("b")
        self.cov_start, self.cov_end, self.cov_w = array("q"), array("q"), array("d")
        self.blocks: dict[tuple, Block] = {}

    def add_interval_cov(self, locus: Locus, value: float) -> None:
        strand = 0 if locus.strand == "+" else 1
        for start, end in locus.intervals:
            self.cov_ref.append(locus.ref_id)
            self.cov_strand.append(strand)
            self.cov_start.append(start)
            self.cov_end.append(end)
            self.cov_w.append(value)

    def add_point(self, locus: Locus, value: float) -> None:
        key = (locus.ref_id, locus.strand, locus.five)
        self.ends[key] = self.ends.get(key, 0.0) + value
        if self.is_5p and locus.cap:
            self.cap[key] = self.cap.get(key, 0.0) + value
            self.cap_positive_reads += value

    def locus_sort_key(self, locus: Locus) -> tuple:
        return (self.refs[locus.ref_id].encode(), locus.five, locus.strand)

    def check_nh(self, records: list[pysam.AlignedSegment], tier: int) -> None:
        """Rule 1: loci found vs NH."""
        nh = next((r.get_tag("NH") for r in records if r.has_tag("NH")), None)
        if nh is not None and nh != tier:
            self.stats["nh_mismatch"] += 1

    def process(self, group: list[pysam.AlignedSegment]) -> None:
        self.stats["reads_total"] += 1
        records = [r for r in group if not r.flag & SKIP_FLAGS]
        if not records:
            fate = "unmapped" if all(r.is_unmapped for r in group) else "discarded_flagged"
            self.stats[fate] += 1
            return

        if any(r.is_paired for r in records):
            if self.compat:
                raise EndCountsError("--nanopare_compat supports single-end reads only (DL-046)")
            pairs = pair_records(
                [r for r in records if not r.is_read2], [r for r in records if r.is_read2]
            )
            loci_all = [pair_locus(m1, m2, self.genome, self.is_5p) for m1, m2 in pairs]
        else:
            loci_all = [record_locus(r, self.genome, self.is_5p) for r in records]

        if self.compat:
            # Rule 4 compat: primary loci first, keyed by start only; colliding records share a
            # locus (#24).
            groups: dict[tuple, list[Locus]] = {}
            for locus in sorted(loci_all, key=lambda x: x.secondary):
                key: tuple = (locus.secondary, locus.aln_start)
                if key in groups:
                    self.stats["loci_key_collisions"] += 1
                groups.setdefault(key, []).append(locus)
            self.check_nh(records, len(groups))
            self.process_compat(records, list(groups.values()))
            return

        # Rule 4: one locus per (reference, strand, 5' end); the first record wins.
        loci: dict[tuple, Locus] = {}
        for locus in loci_all:
            if locus.status == "no_mate1":
                key = ("no_mate1", locus.ref_id, locus.strand, locus.aln_start)
            else:
                key = (locus.ref_id, locus.strand, locus.five)
            loci.setdefault(key, locus)
        tier = len(loci)
        self.check_nh(records, tier)

        # Rules 8 and 11; a locus is counted once, under its first failing reason.
        usable: list[Locus] = []
        for locus in loci.values():
            if locus.status != "no_mate1":
                if locus.n_positions < self.args.min_aligned:
                    locus.status = "too_short"
                elif self.is_5p and locus.untemplated5.strip("G"):
                    locus.status = "untemplated_non_g"
            if locus.status == "ok":
                usable.append(locus)
            else:
                self.stats[f"loci_{locus.status}"] += 1

        if not usable:
            reasons = {locus.status for locus in loci.values()}
            for reason in ("too_short", "untemplated_non_g", "no_mate1"):
                if reason in reasons:
                    self.stats[f"discarded_{reason}"] += 1
                    return

        if self.is_5p and any(locus.has_untemplated3 for locus in usable):
            self.stats["kept_3prime_clip_only"] += 1
        weight = read_weight(records)

        if len(usable) == 1:
            # Rule 14
            locus = usable[0]
            self.stats["unique" if tier == 1 else "unique_after_filter"] += 1
            self.add_point(locus, weight)
            self.add_interval_cov(locus, weight)
            if self.is_5p:
                self.g_runs[g_run_bin(len(locus.untemplated5))] += 1
            return

        if self.args.multimapper_mode == "unique":
            self.stats["multimapped_dropped"] += 1
            return

        # Rule 20: blocks of identical reads with the same locus set (D3).
        usable.sort(key=self.locus_sort_key)
        locus_keys = tuple(self.locus_sort_key(locus) for locus in usable)
        minpos = min(locus.aln_start for locus in loci.values())
        sequence = read_sequence(records)
        block_key = (tier, minpos, sequence, locus_keys)
        block = self.blocks.get(block_key)
        if block is None:
            block = Block(loci=usable, sort_key=(minpos, sequence.encode(), locus_keys))
            self.blocks[block_key] = block
        block.weight += weight
        block.n_reads += 1

    def clip_failure(self, locus: Locus) -> str | None:
        """Compat clip checks over both clips (rule 9, #34)."""
        if not (locus.five_sites or locus.three_sites):
            return None
        if not locus.five_sites:
            return "3prime_clip_only"
        if any(base != "G" for _, base in locus.five_sites + locus.three_sites):
            return "untemplated_non_g"
        return None

    def add_cap_sites(self, locus: Locus, value: float) -> None:
        """Compat cap track: one count per untemplated G at its own coordinate (#8, #34)."""
        length = self.genome.lengths[self.refs[locus.ref_id]]
        for pos, _ in locus.five_sites + locus.three_sites:
            if 0 <= pos < length:
                key = (locus.ref_id, locus.strand, pos)
                self.cap[key] = self.cap.get(key, 0) + value

    def merge_group(self, group: list[Locus], keep: list[Locus] | None = None) -> Locus:
        """One compat locus from records sharing a start: union of positions (#24)."""
        first = group[0]
        members = group if keep is None else keep
        intervals = merge_intervals([iv for locus in members for iv in locus.intervals])
        n = sum(e - s for s, e in intervals)
        five = intervals[0][0] if first.strand == "+" else intervals[-1][1] - 1
        with_sites = [locus for locus in members if locus.five_sites or locus.three_sites]
        last = with_sites[-1] if with_sites else first
        failures = [self.clip_failure(locus) for locus in group] if self.is_5p else []
        # Length check: union of CIGAR positions of all colliding records, before reclaim.
        cigar_ivs = merge_intervals([iv for locus in group for iv in locus.cigar_intervals])
        return Locus(
            ref_id=first.ref_id,
            strand=first.strand,
            five=five,
            intervals=intervals,
            n_positions=n,
            aln_start=first.aln_start,
            status=next((f for f in failures if f), None) or "ok",
            n_cigar=sum(e - s for s, e in cigar_ivs),
            five_sites=last.five_sites,
            three_sites=last.three_sites,
            cigar=min(locus.cigar for locus in group),
            secondary=first.secondary,
        )

    def process_compat(
        self, records: list[pysam.AlignedSegment], groups: list[list[Locus]]
    ) -> None:
        tier = len(groups)
        usable = [g for g in groups if g[0].status != "no_mate1"]
        self.stats["loci_no_mate1"] += len(groups) - len(usable)
        if not usable:
            self.stats["discarded_no_mate1"] += 1
            return
        weight = read_weight(records)

        if tier == 1:
            # Each record is checked on its own; failing records add nothing (#24, #34).
            group = usable[0]
            failures = [self.clip_failure(r) if self.is_5p else None for r in group]
            passing = [r for r, f in zip(group, failures, strict=True) if f is None]
            for r in passing:
                if self.is_5p:
                    # G's recorded before the length check (#33), unweighted int (#38).
                    self.add_cap_sites(r, 1)
            locus = self.merge_group(group, passing) if passing else None
            too_short = locus is None or locus.n_positions <= self.args.min_aligned
            if too_short:
                failure = next((f for f in failures if f), None)
                reason = failure if locus is None and failure else "too_short"
                if reason != "3prime_clip_only":
                    self.stats[f"loci_{reason}"] += 1
                self.stats[f"discarded_{reason}"] += 1
                return
            assert locus is not None
            self.stats["unique"] += 1
            key = (locus.ref_id, locus.strand, locus.five)
            self.ends[key] = self.ends.get(key, 0.0) + weight
            self.add_interval_cov(locus, weight)
            if self.is_5p:
                self.g_runs[g_run_bin(len(locus.five_sites))] += 1
                if locus.five_sites:
                    self.cap_positive_reads += weight
            return

        # Block of identical reads; it keeps the last read's loci and weight (G1, #36).
        merged = [self.merge_group(g) for g in usable]
        minpos = min(locus.aln_start for locus in merged)
        if self.is_5p:
            ident = read_sequence(records)
        else:
            ident = next(locus.cigar for locus in merged if locus.aln_start == minpos)
        usable = merged
        block_key = (tier, minpos, ident)
        block = self.blocks.get(block_key)
        if block is None:
            block = Block(loci=usable, sort_key=(minpos,))
            self.blocks[block_key] = block
        block.loci = usable
        block.weight = weight
        block.n_reads += 1

    def share_compat(self, block: Block, cov: dict, final: bool) -> bool:
        """Rule 20 compat for one block; False if it stays for the next pass."""
        keep = [i for i, locus in enumerate(block.loci) if locus.n_cigar >= self.args.min_aligned]
        if not keep:
            self.stats["loci_too_short"] += len(block.loci)
            self.stats["discarded_too_short"] += block.n_reads
            return True
        loci = [block.loci[i] for i in keep]
        idx = [block.shares_idx[i] for i in keep]
        weights = [
            float(cov[(locus.ref_id, locus.strand)][i].sum())
            for locus, i in zip(loci, idx, strict=True)
        ]
        if not self.is_5p and sum(weights) == 0:
            flip = {"+": "-", "-": "+"}
            weights = [
                float(cov[(locus.ref_id, flip[locus.strand])][i].sum())
                for locus, i in zip(loci, idx, strict=True)
            ]
        naive = sum(weights) == 0
        if naive:
            if not final:
                return False
            weights = [1.0] * len(loci)
        total = float(sum(weights))
        shares = [w / total for w in weights]
        if self.is_5p:
            for i, locus in enumerate(loci):
                if locus.status != "ok":
                    shares[i] = 0.0
                    if locus.status == "untemplated_non_g":
                        self.stats["loci_untemplated_non_g"] += 1
        if sum(shares) == 0:
            self.stats["multimapped_lost_compat"] += block.n_reads
            return True
        norm = sum(shares)
        shares = [s / norm for s in shares]
        for locus, i, share in zip(loci, idx, shares, strict=True):
            if share == 0:
                continue
            value = float(block.n_reads) * share * block.weight
            key = (locus.ref_id, locus.strand, locus.five)
            self.ends[key] = self.ends.get(key, 0.0) + value
            self.add_interval_cov(locus, value)
            cov[(locus.ref_id, locus.strand)][i] += value
            if self.is_5p and locus.five_sites:
                self.add_cap_sites(locus, value)
                self.cap_positive_reads += value
        self.stats["multimapped_naive" if naive else "multimapped_rescued"] += block.n_reads
        return True

    def run_multimappers_compat(self) -> None:
        if not self.blocks:
            return
        by_tier: dict[int, list[Block]] = {}
        for key, block in self.blocks.items():
            by_tier.setdefault(key[0], []).append(block)
        cov = self.rescue_coverage(by_tier)
        for tier in sorted(by_tier):
            # Blocks per leftmost start, first-seen order within a start.
            by_start: dict[int, list[Block]] = {}
            for block in by_tier[tier]:
                by_start.setdefault(block.sort_key[0], []).append(block)
            ascending = True
            cleared = True
            while cleared:
                cleared = False
                for start in sorted(by_start, reverse=not ascending):
                    left = [b for b in by_start[start] if not self.share_compat(b, cov, False)]
                    if left:
                        by_start[start] = left
                    else:
                        del by_start[start]
                        cleared = True
                ascending = not ascending
            for start in sorted(by_start):
                for block in by_start[start]:
                    self.share_compat(block, cov, True)

    def place(self, block: Block, shares: list[float], cov: dict | None) -> None:
        """Rule 19; `cov` is updated at once so later blocks see this one."""
        for i, (locus, share) in enumerate(zip(block.loci, shares, strict=True)):
            value = block.weight * share
            self.add_point(locus, value)
            self.add_interval_cov(locus, value)
            if cov is not None:
                cov[(locus.ref_id, locus.strand)][block.shares_idx[i]] += value

    def run_multimappers(self) -> None:
        if self.compat:
            self.run_multimappers_compat()
            return
        if not self.blocks:
            return
        by_tier: dict[int, list[Block]] = {}
        for key, block in self.blocks.items():
            by_tier.setdefault(key[0], []).append(block)
        for tier_blocks in by_tier.values():
            tier_blocks.sort(key=lambda b: b.sort_key)

        if self.args.multimapper_mode == "uniform":
            for tier in sorted(by_tier):
                for block in by_tier[tier]:
                    k = len(block.loci)
                    self.place(block, [1.0 / k] * k, None)
                    self.stats["multimapped_uniform"] += block.n_reads
            return

        cov = self.rescue_coverage(by_tier)
        for tier in sorted(by_tier):
            remaining = by_tier[tier]
            ascending = True
            while remaining:
                kept: list[Block] = []
                order = remaining if ascending else remaining[::-1]
                for block in order:
                    weights = self.rescue_weights(block, cov)
                    if weights is None:
                        kept.append(block)
                        continue
                    total = sum(weights)
                    self.place(block, [w / total for w in weights], cov)
                    self.stats["multimapped_rescued"] += block.n_reads
                if len(kept) == len(remaining):
                    break
                remaining = kept if ascending else kept[::-1]
                ascending = not ascending
            for block in remaining:
                k = len(block.loci)
                self.place(block, [1.0 / k] * k, cov)
                self.stats["multimapped_naive"] += block.n_reads

    def rescue_weights(self, block: Block, cov: dict) -> list[float] | None:
        """Rule 18; None when every locus has weight 0."""
        weights = [
            float(cov[(locus.ref_id, locus.strand)][idx].sum())
            for locus, idx in zip(block.loci, block.shares_idx, strict=True)
        ]
        if sum(weights) > 0:
            return weights
        if not self.is_5p:
            flip = {"+": "-", "-": "+"}
            weights = [
                float(cov[(locus.ref_id, flip[locus.strand])][idx].sum())
                for locus, idx in zip(block.loci, block.shares_idx, strict=True)
            ]
            if sum(weights) > 0:
                return weights
        return None

    def rescue_coverage(self, by_tier: dict[int, list[Block]]) -> dict:
        """Unique-read coverage at the positions multimapper loci cover, one array per strand."""
        positions: dict[int, list[np.ndarray]] = {}
        for blocks in by_tier.values():
            for block in blocks:
                for locus in block.loci:
                    for s, e in locus.intervals:
                        positions.setdefault(locus.ref_id, []).append(np.arange(s, e))
        grid = {ref: np.unique(np.concatenate(parts)) for ref, parts in positions.items()}
        cov = {}
        for ref, p in grid.items():
            for strand in "+-":
                cov[(ref, strand)] = np.zeros(len(p), dtype=np.float64)

        refs = np.frombuffer(self.cov_ref, dtype=np.int32)
        strands = np.frombuffer(self.cov_strand, dtype=np.int8)
        starts = np.frombuffer(self.cov_start, dtype=np.int64)
        ends = np.frombuffer(self.cov_end, dtype=np.int64)
        weights = np.frombuffer(self.cov_w, dtype=np.float64)
        for ref, p in grid.items():
            for strand_code, strand in ((0, "+"), (1, "-")):
                mask = (refs == ref) & (strands == strand_code)
                i0 = np.searchsorted(p, starts[mask], "left")
                i1 = np.searchsorted(p, ends[mask], "left")
                lengths = i1 - i0
                total = int(lengths.sum())
                if total == 0:
                    continue
                offsets = np.repeat(np.cumsum(lengths) - lengths - i0, lengths)
                idx = np.arange(total) - offsets
                # add.at is unbuffered: adds happen in BAM order (rule 24).
                np.add.at(cov[(ref, strand)], idx, np.repeat(weights[mask], lengths))

        for blocks in by_tier.values():
            for block in blocks:
                block.shares_idx = [
                    np.searchsorted(
                        grid[locus.ref_id],
                        np.concatenate([np.arange(s, e) for s, e in locus.intervals]),
                    )
                    for locus in block.loci
                ]
        return cov

    def write_compat(self) -> None:
        """Rule 23 compat and DL-044: the original's bedgraph format."""
        prefix = self.args.prefix
        order = sorted(range(len(self.refs)), key=lambda i: self.refs[i].encode())
        tracks = [("ends", self.ends, False)] + ([("cap", self.cap, True)] if self.is_5p else [])
        for name, track, cap in tracks:
            for strand, label in (("+", "plus"), ("-", "minus")):
                by_ref: dict[int, dict[int, float]] = {}
                for (ref_id, s, pos), value in track.items():
                    if s == strand:
                        by_ref.setdefault(ref_id, {})[pos] = value
                with open(f"{prefix}.{name}.{label}.bedgraph", "w") as out:
                    for ref_id in order:
                        if ref_id in by_ref:
                            write_compat_runs(out, self.refs[ref_id], by_ref[ref_id], cap)
        if self.args.write_coverage:
            # One reference and strand in memory at a time.
            index = self.coverage_index()
            if self.is_5p:
                for code, label in ((0, "plus"), (1, "minus")):
                    with open(f"{prefix}.cov.{label}.bedgraph", "w") as out:
                        for ref_id in order:
                            if (ref_id, code) in index:
                                values = self.coverage_array(ref_id, index[(ref_id, code)])
                                write_compat_runs(
                                    out, self.refs[ref_id], nonzero_dict(values), False
                                )
            else:
                # #39: strand values already rounded to 2 dp, then summed.
                with open(f"{prefix}.cov.bedgraph", "w") as out:
                    for ref_id in order:
                        summed = np.zeros(self.genome.lengths[self.refs[ref_id]], dtype=np.float64)
                        for code in (0, 1):
                            if (ref_id, code) not in index:
                                continue
                            values = self.coverage_array(ref_id, index[(ref_id, code)])
                            nz = np.flatnonzero(values)
                            summed[nz] += [round(v, 2) for v in values[nz].tolist()]
                        if summed.any():
                            write_compat_runs(out, self.refs[ref_id], nonzero_dict(summed), False)

    def coverage_index(self) -> dict[tuple[int, int], np.ndarray]:
        """Coverage interval indices per (reference, strand code), in read order."""
        refs = np.frombuffer(self.cov_ref, dtype=np.int32)
        strands = np.frombuffer(self.cov_strand, dtype=np.int8)
        if refs.size == 0:
            return {}
        order = np.lexsort((strands, refs))  # stable: read order kept within a group
        keys = refs[order].astype(np.int64) * 2 + strands[order]
        groups = np.split(order, np.flatnonzero(np.diff(keys)) + 1)
        return {(int(refs[g[0]]), int(strands[g[0]])): g for g in groups}

    def coverage_array(self, ref_id: int, idx: np.ndarray) -> np.ndarray:
        """Per-nt coverage of one reference and strand."""
        starts = np.frombuffer(self.cov_start, dtype=np.int64)[idx]
        ends = np.frombuffer(self.cov_end, dtype=np.int64)[idx]
        weights = np.frombuffer(self.cov_w, dtype=np.float64)[idx]
        values = np.zeros(self.genome.lengths[self.refs[ref_id]], dtype=np.float64)
        # Per-position adds in read order, as the original sums floats (no cumsum error);
        # chunked so the expanded positions stay small.
        for i in range(0, len(idx), COVERAGE_CHUNK):
            s, e = starts[i : i + COVERAGE_CHUNK], ends[i : i + COVERAGE_CHUNK]
            lengths = e - s
            offsets = np.repeat(np.cumsum(lengths) - lengths - s, lengths)
            positions = np.arange(int(lengths.sum())) - offsets
            chunk_weights = np.repeat(weights[i : i + COVERAGE_CHUNK], lengths)
            # Past the reference end only after a cross-chromosome collision (#24): skipped.
            inside = positions < values.size
            np.add.at(values, positions[inside], chunk_weights[inside])
        return values

    def write(self) -> None:
        if self.compat:
            self.write_compat()
            return
        prefix, decimals = self.args.prefix, self.args.decimals
        order = sorted(range(len(self.refs)), key=lambda i: self.refs[i].encode())
        rank = {ref_id: n for n, ref_id in enumerate(order)}
        tracks = [("ends", self.ends)] + ([("cap", self.cap)] if self.is_5p else [])
        for name, track in tracks:
            for strand, label in (("+", "plus"), ("-", "minus")):
                keys = sorted(
                    (k for k in track if k[1] == strand), key=lambda k: (rank[k[0]], k[2])
                )
                with open(f"{prefix}.{name}.{label}.bedgraph", "w") as out:
                    for ref_id, _, pos in keys:
                        value = round(track[(ref_id, strand, pos)], decimals)
                        if value != 0:
                            text = format_value(value, decimals)
                            out.write(f"{self.refs[ref_id]}\t{pos}\t{pos + 1}\t{text}\n")
        if self.args.write_coverage:
            self.write_coverage(order)

    def write_coverage(self, order: list[int]) -> None:
        refs = np.frombuffer(self.cov_ref, dtype=np.int32)
        strands = np.frombuffer(self.cov_strand, dtype=np.int8)
        starts = np.frombuffer(self.cov_start, dtype=np.int64)
        ends = np.frombuffer(self.cov_end, dtype=np.int64)
        weights = np.frombuffer(self.cov_w, dtype=np.float64)
        prefix, decimals = self.args.prefix, self.args.decimals
        if self.is_5p:
            outputs = {0: [0], 1: [1]}
            names = {0: f"{prefix}.cov.plus.bedgraph", 1: f"{prefix}.cov.minus.bedgraph"}
        else:
            outputs = {0: [0, 1]}
            names = {0: f"{prefix}.cov.bedgraph"}
        for out_id, strand_codes in outputs.items():
            with open(names[out_id], "w") as out:
                for ref_id in order:
                    mask = (refs == ref_id) & np.isin(strands, strand_codes)
                    if not mask.any():
                        continue
                    length = self.genome.lengths[self.refs[ref_id]]
                    diff = np.zeros(length + 1, dtype=np.float64)
                    np.add.at(diff, starts[mask], weights[mask])
                    np.add.at(diff, ends[mask], -weights[mask])
                    values = np.round(np.cumsum(diff[:-1]), decimals)
                    write_runs(out, self.refs[ref_id], values, decimals)

    def stats_json(self) -> dict:
        stats: dict = {
            "schema_version": SCHEMA_VERSION,
            "parameters": {
                "library": self.args.library,
                "min_aligned": self.args.min_aligned,
                "multimapper_mode": self.args.multimapper_mode,
                "decimals": self.args.decimals,
                "nanopare_compat": self.args.nanopare_compat,
            },
            "reads_total": self.stats["reads_total"],
            "read_fates": {k: self.stats[k] for k in READ_FATES},
            "loci": {k: self.stats[k] for k in LOCUS_COUNTERS},
            "kept_3prime_clip_only": self.stats["kept_3prime_clip_only"],
            "nh_mismatch": self.stats["nh_mismatch"],
        }
        if self.is_5p:
            stats["cap_positive_reads"] = round(self.cap_positive_reads, self.args.decimals)
            stats["untemplated_g_run_length"] = self.g_runs
        return stats


def write_runs(out, chrom: str, values: np.ndarray, decimals: int) -> None:
    """Coverage rows: runs of equal non-zero value merged."""
    if values.size == 0:
        return
    change = np.flatnonzero(np.diff(values)) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [values.size]))
    for s, e in zip(starts.tolist(), ends.tolist(), strict=True):
        value = float(values[s])
        if value != 0:
            out.write(f"{chrom}\t{s}\t{e}\t{format_value(value, decimals)}\n")


def nonzero_dict(values: np.ndarray) -> dict[int, float]:
    idx = np.flatnonzero(values)
    return dict(zip(idx.tolist(), values[idx].tolist(), strict=True))


def write_compat_runs(out, chrom: str, values: dict[int, float], cap: bool) -> None:
    """Compat rows: rounded to 2 dp, then equal adjacent values merged.

    Cap rows are compared as text and written even when `0.0` (#38); a chromosome/strand whose
    only value is at position 0 is not written (#37).
    """
    if not values or list(values) == [0]:
        return
    run: list | None = None
    for pos in sorted(values):
        text = str(round(values[pos], 2))
        keep = cap or float(text) != 0
        if run is not None and keep and pos == run[1] and text == run[2]:
            run[1] = pos + 1
            continue
        if run is not None:
            out.write(f"{chrom}\t{run[0]}\t{run[1]}\t{run[2]}\n")
        run = [pos, pos + 1, text] if keep else None
    if run is not None:
        out.write(f"{chrom}\t{run[0]}\t{run[1]}\t{run[2]}\n")


def open_inputs(args: argparse.Namespace) -> tuple[pysam.AlignmentFile, Genome, list[str]]:
    bam = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
    if bam.header.to_dict().get("HD", {}).get("SO") == "coordinate":
        raise EndCountsError("BAM is coordinate-sorted; reads must be grouped by name (rule 1)")
    refs = list(bam.references)
    lengths = dict(zip(refs, bam.lengths, strict=True))
    fasta = None
    if args.fasta is not None:
        fasta = pysam.FastaFile(args.fasta)
        fasta_lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
        for ref, length in lengths.items():
            if fasta_lengths.get(ref) != length:
                raise EndCountsError(f"reference {ref} missing from FASTA or length differs")
    return bam, Genome(fasta, lengths), refs


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    if args.nanopare_compat and args.write_coverage:
        print(
            "degrseq end-counts: warning: --write_coverage with --nanopare_compat keeps every "
            "aligned interval in memory (about 30 bytes each) plus up to two arrays the length "
            "of the longest chromosome (8 bytes per base each); request memory accordingly",
            file=sys.stderr,
        )
    try:
        bam, genome, refs = open_inputs(args)
        counter = Counter(args, refs, genome)
        with bam:
            for group in read_groups(bam):
                counter.process(group)
        counter.run_multimappers()
    except EndCountsError as error:
        sys.exit(f"degrseq end-counts: error: {error}")
    counter.write()
    stats_path = Path(f"{args.prefix}.end-counts.stats.json")
    stats_path.write_text(json.dumps(counter.stats_json(), indent=2) + "\n")
