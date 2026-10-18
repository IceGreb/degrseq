"""Unit tests for specs/end-counts.md, default mode. Each test names the rule it checks."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pysam
import pytest
from end_counts_helpers import (
    CHROM_LENGTH,
    GENOME,
    build_inputs,
    fates,
    run,
    sam,
)

from degrseq import end_counts

G1 = GENOME["chr1"]


def first_pos(predicate, start: int = 30) -> int:
    return next(p for p in range(start, CHROM_LENGTH - 40) if predicate(p))


def unique(name: str, pos: int, chrom: str = "chr1", flag: int = 0, cigar: str = "20M") -> str:
    return sam(name, flag, chrom, pos, cigar, tags="NH:i:1")


def multi(name: str, loci: list[tuple], seq: str | None = None, tags: str = "") -> list[str]:
    """One record per locus (chrom, pos[, flag[, cigar]]); the first is primary."""
    records = []
    for i, locus in enumerate(loci):
        chrom, pos = locus[0], locus[1]
        flag = (locus[2] if len(locus) > 2 else 0) | (256 if i else 0)
        cigar = locus[3] if len(locus) > 3 else "20M"
        extra = f"NH:i:{len(loci)}\tHI:i:{i + 1}" + (f"\t{tags}" if tags else "")
        records.append(sam(name, flag, chrom, pos, cigar, seq, extra))
    return records


# Rule 1
def test_coordinate_sorted_bam_is_rejected(tmp_path: Path) -> None:
    bam, fasta = build_inputs(tmp_path, [unique("r1", 50)], sort_order="coordinate")
    with pytest.raises(SystemExit, match="coordinate-sorted"):
        end_counts.main(
            ["--bam", str(bam), "--fasta", str(fasta), "--library", "5P", "--prefix", "x"]
        )


def test_nh_mismatch_is_counted_and_read_kept(tmp_path: Path) -> None:
    result = run(tmp_path, [sam("r1", 0, "chr1", 50, "20M", tags="NH:i:2")])
    assert result["stats"]["nh_mismatch"] == 1
    assert fates(result) == {"unique": 1}


# Rule 2
def test_flag_filtering(tmp_path: Path) -> None:
    records = [
        sam("unmapped", 4, "*", -1, "*", "A" * 20),
        sam("dup", 1024, "chr1", 50, "20M"),
        sam("qcfail", 512, "chr1", 50, "20M"),
        sam("supp", 0, "chr1", 50, "20M", tags="NH:i:1"),
        sam("supp", 2048, "chr1", 150, "20M"),
        *multi("secondary", [("chr1", 200), ("chr2", 200)]),
    ]
    result = run(tmp_path, records, extra=["--multimapper_mode", "uniform"])
    assert fates(result) == {
        "unmapped": 1,
        "discarded_flagged": 2,
        "unique": 1,
        "multimapped_uniform": 1,
    }
    assert result[("ends", "+")] == {
        ("chr1", 50): "1",
        ("chr1", 200): "0.5",
        ("chr2", 200): "0.5",
    }


# Rule 4
def test_same_position_on_two_chromosomes_is_two_loci(tmp_path: Path) -> None:
    result = run(tmp_path, multi("r1", [("chr1", 60), ("chr2", 60)]))
    assert fates(result) == {"multimapped_naive": 1}
    assert result[("ends", "+")] == {("chr1", 60): "0.5", ("chr2", 60): "0.5"}


def test_records_with_one_locus_key_count_once(tmp_path: Path) -> None:
    # Same 5' end, different 3' extent: one locus, so the read is unique.
    result = run(tmp_path, multi("r1", [("chr1", 60, 0, "20M"), ("chr1", 60, 0, "18M")]))
    assert fates(result) == {"unique": 1}
    assert result[("ends", "+")] == {("chr1", 60): "1"}


# Rules 6, 7, 10
def test_reclaim_then_untemplated_g_is_cap(tmp_path: Path) -> None:
    p = 80
    seq = "G" + G1[p - 2 : p] + G1[p : p + 20]
    result = run(tmp_path, [sam("r1", 0, "chr1", p, "3S20M", seq, "NH:i:1")])
    assert result[("ends", "+")] == {("chr1", p - 2): "1"}
    assert result[("cap", "+")] == {("chr1", p - 2): "1"}
    assert result["stats"]["untemplated_g_run_length"]["1"] == 1


def test_reclaim_stops_at_first_mismatch(tmp_path: Path) -> None:
    # Clip x G y with y = genome[p-1] and x = genome[p-3]: only y is reclaimed, so x stays
    # untemplated and the locus is unusable even though x matches the genome.
    p = 80
    seq = G1[p - 3] + "G" + G1[p - 1] + G1[p : p + 20]
    result = run(tmp_path, [sam("r1", 0, "chr1", p, "3S20M", seq, "NH:i:1")])
    assert fates(result) == {"discarded_untemplated_non_g": 1}
    assert result["stats"]["loci"]["loci_untemplated_non_g"] == 1


def test_reclaim_stops_at_reference_start(tmp_path: Path) -> None:
    seq = "GG" + G1[0:2] + G1[2:22]
    result = run(tmp_path, [sam("r1", 0, "chr1", 2, "4S20M", seq, "NH:i:1")])
    assert result[("cap", "+")] == {("chr1", 0): "1"}
    assert result["stats"]["untemplated_g_run_length"]["2"] == 1


def test_reclaim_stops_at_reference_end_minus_strand(tmp_path: Path) -> None:
    p = CHROM_LENGTH - 22
    seq = G1[p : p + 20] + G1[p + 20 : p + 22] + "C"  # 3' clip in reference orientation
    result = run(tmp_path, [sam("r1", 16, "chr1", p, "20M3S", seq, "NH:i:1")])
    assert result[("ends", "-")] == {("chr1", CHROM_LENGTH - 1): "1"}
    assert result[("cap", "-")] == {("chr1", CHROM_LENGTH - 1): "1"}


def test_genome_n_is_never_reclaimed(tmp_path: Path) -> None:
    # chr2:100 is N; a read N there must stay untemplated, so the locus has a non-G 5' base.
    g2 = GENOME["chr2"]
    seq = "GN" + g2[101:121]
    result = run(tmp_path, [sam("r1", 0, "chr2", 101, "2S20M", seq, "NH:i:1")])
    assert fates(result) == {"discarded_untemplated_non_g": 1}


def test_five_prime_end_minus_strand_no_clip(tmp_path: Path) -> None:
    result = run(tmp_path, [unique("r1", 50, flag=16)])
    assert result[("ends", "-")] == {("chr1", 69): "1"}
    assert result[("ends", "+")] == {}


def test_minus_strand_ggg_cap(tmp_path: Path) -> None:
    p = first_pos(lambda q: G1[q + 20] != "C")
    seq = G1[p : p + 20] + "CCC"
    result = run(tmp_path, [sam("r1", 16, "chr1", p, "20M3S", seq, "NH:i:1")])
    assert result[("ends", "-")] == {("chr1", p + 19): "1"}
    assert result[("cap", "-")] == {("chr1", p + 19): "1"}
    assert result["stats"]["untemplated_g_run_length"]["3"] == 1


# Rule 8
def test_untemplated_ga_makes_locus_unusable(tmp_path: Path) -> None:
    p = first_pos(lambda q: G1[q - 1] != "A")
    seq = "GA" + G1[p : p + 20]
    result = run(tmp_path, [sam("r1", 0, "chr1", p, "2S20M", seq, "NH:i:1")])
    assert fates(result) == {"discarded_untemplated_non_g": 1}
    assert result[("ends", "+")] == {}


# Rule 9
def test_three_prime_only_clip_is_kept(tmp_path: Path) -> None:
    p = 60
    seq = G1[p : p + 20] + "GG"
    result = run(tmp_path, [sam("r1", 0, "chr1", p, "20M2S", seq, "NH:i:1")])
    assert fates(result) == {"unique": 1}
    assert result["stats"]["kept_3prime_clip_only"] == 1
    assert result[("cap", "+")] == {}


# Rule 11
def test_min_aligned_boundary(tmp_path: Path) -> None:
    records = [unique("len16", 50, cigar="16M"), unique("len15", 100, cigar="15M")]
    result = run(tmp_path, records)
    assert fates(result) == {"unique": 1, "discarded_too_short": 1}
    assert result[("ends", "+")] == {("chr1", 50): "1"}


def test_reclaimed_bases_count_towards_min_aligned_for_multimappers(tmp_path: Path) -> None:
    # 15 aligned + 1 reclaimed = 16 at both loci.
    p = 60
    q = first_pos(lambda x: G1[x - 1] == G1[p - 1], start=120)
    loci = [("chr1", p, 0, "1S15M"), ("chr1", q, 0, "1S15M")]
    records = multi("r1", loci, seq=G1[p - 1] + G1[p : p + 15])
    assert fates(run(tmp_path / "a", records)) == {"multimapped_naive": 1}
    stricter = run(tmp_path / "b", records, extra=["--min_aligned", "17"])
    assert fates(stricter) == {"discarded_too_short": 1}


def test_tier_counts_loci_before_length_filter(tmp_path: Path) -> None:
    records = multi("r1", [("chr1", 60, 0, "20M"), ("chr1", 160, 0, "10M")])
    result = run(tmp_path, records)
    assert fates(result) == {"unique_after_filter": 1}
    assert result["stats"]["loci"]["loci_too_short"] == 1


# Rule 13
def test_xw_weight_applies_to_all_tracks(tmp_path: Path) -> None:
    seq = "GGG" + G1[50:70]
    records = [sam("r1", 0, "chr1", 50, "3S20M", seq, "NH:i:1\tXW:f:0.5")]
    result = run(tmp_path, records, extra=["--write_coverage"])
    assert result[("ends", "+")] == {("chr1", 50): "0.5"}
    assert result[("cap", "+")] == {("chr1", 50): "0.5"}
    assert Path(f"{result['prefix']}.cov.plus.bedgraph").read_text() == "chr1\t50\t70\t0.5\n"
    assert result["stats"]["cap_positive_reads"] == 0.5


def test_conflicting_xw_is_an_error(tmp_path: Path) -> None:
    records = multi("r1", [("chr1", 60), ("chr1", 160)])
    records[0] += "\tXW:f:0.5"
    records[1] += "\tXW:f:1"
    with pytest.raises(SystemExit, match="different XW"):
        run(tmp_path, records)


# Rules 18, 19
def test_rescue_shares_follow_coverage(tmp_path: Path) -> None:
    records = [*[unique(f"u{i}", 60) for i in range(3)], unique("u3", 160)]
    records += multi("m", [("chr1", 60), ("chr1", 160)])
    result = run(tmp_path, records)
    assert result[("ends", "+")] == {("chr1", 60): "3.75", ("chr1", 160): "1.25"}
    assert fates(result) == {"unique": 4, "multimapped_rescued": 1}


def test_no_coverage_gives_equal_split(tmp_path: Path) -> None:
    result = run(tmp_path, multi("m", [("chr1", 60), ("chr1", 160)]))
    assert result[("ends", "+")] == {("chr1", 60): "0.5", ("chr1", 160): "0.5"}
    assert fates(result) == {"multimapped_naive": 1}


def test_body_falls_back_to_opposite_strand(tmp_path: Path) -> None:
    records = [unique("u", 10, flag=16), *multi("m", [("chr1", 10), ("chr1", 200)])]
    result = run(tmp_path, records, library="BODY")
    assert result[("ends", "+")] == {("chr1", 10): "1"}
    assert fates(result) == {"unique": 1, "multimapped_rescued": 1}


def test_5p_has_no_opposite_strand_fallback(tmp_path: Path) -> None:
    records = [unique("u", 10, flag=16), *multi("m", [("chr1", 10), ("chr1", 200)])]
    result = run(tmp_path, records)
    assert result[("ends", "+")] == {("chr1", 10): "0.5", ("chr1", 200): "0.5"}


# Rule 20
def test_ascending_first_pass_with_immediate_update(tmp_path: Path) -> None:
    records = [unique("u1", 100), unique("u2", 100), unique("u3", 200)]
    records += multi("m3", [("chr1", 200), ("chr2", 250)], seq="C" * 20)  # minpos 200
    records += multi("m1", [("chr1", 100), ("chr1", 200)], seq="A" * 20)  # minpos 100
    result = run(tmp_path, records)
    # m1 first: 2/3 : 1/3; then m3 sees chr1:200 covered and chr2:250 empty.
    assert result[("ends", "+")] == {("chr1", 100): "2.6667", ("chr1", 200): "2.3333"}


def test_zero_coverage_block_deferred_to_next_pass(tmp_path: Path) -> None:
    records = [unique("u", 90)]
    records += multi("d1", [("chr1", 10), ("chr1", 115)], seq="A" * 20)  # minpos 10, no cov
    records += multi("d2", [("chr1", 100), ("chr1", 250)], seq="C" * 20)  # overlaps u
    result = run(tmp_path, records)
    ends = result[("ends", "+")]
    assert ends[("chr1", 115)] == "1"
    assert ("chr1", 10) not in ends
    assert fates(result) == {"unique": 1, "multimapped_rescued": 2}


def test_identical_reads_form_one_block(tmp_path: Path) -> None:
    # Unequal footprints (25 vs 20 nt), so sharing twice in sequence would differ from once.
    records = [unique("u1", 0, cigar="25M"), unique("u2", 100, cigar="25M")]
    loci = [("chr1", 0, 0, "25M"), ("chr1", 100, 0, "5S20M")]
    for name in ("m1", "m2"):
        records += multi(name, loci, seq="A" * 25)
    result = run(tmp_path, records, library="BODY")
    assert result[("ends", "+")] == {("chr1", 0): "2.1111", ("chr1", 100): "1.8889"}


def test_identical_reads_with_different_locus_sets_stay_separate(tmp_path: Path) -> None:
    # D3 (DL-033): b1 keeps its own locus chr1:150; the original would move it to chr1:250.
    records = multi("b1", [("chr1", 50), ("chr1", 150)], seq="A" * 20)
    records += multi("b2", [("chr1", 50), ("chr1", 250)], seq="A" * 20)
    result = run(tmp_path, records)
    ends = result[("ends", "+")]
    assert ("chr1", 150) in ends
    assert ("chr1", 250) in ends


def tie_records(order: list[str]) -> list[str]:
    blocks = {
        "t1": multi("t1", [("chr1", 10), ("chr2", 50)], seq="C" * 20),
        "t2": multi("t2", [("chr2", 10), ("chr2", 50)], seq="A" * 20),
    }
    records = [unique("u1", 10), unique("u2", 10, chrom="chr2"), unique("u3", 50, chrom="chr2")]
    for name in order:
        records += blocks[name]
    return records


def test_tie_order_is_by_sequence_not_bam_order(tmp_path: Path) -> None:
    a = run(tmp_path / "a", tie_records(["t1", "t2"]))
    b = run(tmp_path / "b", tie_records(["t2", "t1"]))
    assert a[("ends", "+")] == b[("ends", "+")]
    # t2 (A...) sorts before t1 (C...): t2 0.5/0.5, then t1 0.4 chr1:10 / 0.6 chr2:50.
    assert a[("ends", "+")][("chr1", 10)] == "1.4"
    assert a[("ends", "+")][("chr2", 10)] == "1.5"
    assert a[("ends", "+")][("chr2", 50)] == "2.1"


# Rule 16
def test_uniform_mode_ignores_coverage(tmp_path: Path) -> None:
    records = [unique(f"u{i}", 60) for i in range(3)]
    records += multi("m", [("chr1", 60), ("chr1", 160)])
    result = run(tmp_path, records, extra=["--multimapper_mode", "uniform"])
    assert result[("ends", "+")] == {("chr1", 60): "3.5", ("chr1", 160): "0.5"}
    assert fates(result) == {"unique": 3, "multimapped_uniform": 1}


def test_unique_mode_drops_multimappers(tmp_path: Path) -> None:
    records = [unique("u", 60), *multi("m", [("chr1", 60), ("chr1", 160)])]
    result = run(tmp_path, records, extra=["--multimapper_mode", "unique"])
    assert result[("ends", "+")] == {("chr1", 60): "1"}
    assert fates(result) == {"unique": 1, "multimapped_dropped": 1}


# Rule 17
def test_paired_end_five_prime_end_from_mate1(tmp_path: Path) -> None:
    records = [
        sam("p", 83, "chr1", 150, "20M", tags="NH:i:1\tHI:i:1"),
        sam("p", 163, "chr1", 100, "20M", tags="NH:i:1\tHI:i:1"),
    ]
    result = run(tmp_path, records, extra=["--write_coverage"])
    assert result[("ends", "-")] == {("chr1", 169): "1"}
    cov = Path(f"{result['prefix']}.cov.minus.bedgraph").read_text()
    assert cov == "chr1\t100\t120\t1\nchr1\t150\t170\t1\n"


def test_paired_end_locus_without_mate1_is_unusable(tmp_path: Path) -> None:
    records = [
        sam("p", 137, "chr1", 100, "20M", tags="NH:i:1\tHI:i:1"),
        sam("p", 69, "*", -1, "*", "A" * 20),
    ]
    result = run(tmp_path, records)
    assert fates(result) == {"discarded_no_mate1": 1}


# Rule 22
def test_ggg_read_counts_one_cap_at_five_prime_end(tmp_path: Path) -> None:
    seq = "GGG" + G1[50:70]
    result = run(tmp_path, [sam("r1", 0, "chr1", 50, "3S20M", seq, "NH:i:1")])
    assert result[("cap", "+")] == {("chr1", 50): "1"}


# Rule 23
def test_rounding_and_format(tmp_path: Path) -> None:
    records = multi("m", [("chr1", 60), ("chr1", 160), ("chr2", 60)])
    records += [sam("w", 0, "chr2", 200, "20M", tags="NH:i:1\tXW:f:0.0033333")]
    four = run(tmp_path / "four", records)
    assert four[("ends", "+")][("chr1", 60)] == "0.3333"
    assert four[("ends", "+")][("chr2", 200)] == "0.0033"
    two = run(tmp_path / "two", records, extra=["--decimals", "2"])
    assert two[("ends", "+")][("chr1", 60)] == "0.33"
    assert ("chr2", 200) not in two[("ends", "+")]


# Rule 24
def test_output_identical_under_different_hash_seeds(tmp_path: Path) -> None:
    records = [*tie_records(["t1", "t2"]), *multi("m", [("chr1", 60), ("chr1", 160)])]
    bam, fasta = build_inputs(tmp_path, records)
    outputs = []
    for seed in ("0", "12345"):
        prefix = tmp_path / f"seed{seed}"
        env = {**os.environ, "PYTHONHASHSEED": seed}
        subprocess.run(
            [sys.executable, "-m", "degrseq", "end-counts", "--bam", str(bam), "--fasta",
             str(fasta), "--library", "5P", "--prefix", str(prefix)],
            check=True,
            env=env,
        )  # fmt: skip
        outputs.append(
            [Path(f"{prefix}.{name}").read_bytes() for name in (
                "ends.plus.bedgraph", "ends.minus.bedgraph", "cap.plus.bedgraph",
                "cap.minus.bedgraph", "end-counts.stats.json")]
        )  # fmt: skip
    assert outputs[0] == outputs[1]


# Rule 25 and BODY coverage
def test_body_coverage_sums_strands_and_merges_runs(tmp_path: Path) -> None:
    records = [unique("a", 10), unique("b", 10, flag=16), unique("c", 20)]
    result = run(tmp_path, records, library="BODY", extra=["--write_coverage"])
    cov = Path(f"{result['prefix']}.cov.bedgraph").read_text()
    assert cov == "chr1\t10\t20\t2\nchr1\t20\t30\t3\nchr1\t30\t40\t1\n"


def test_no_coverage_files_by_default(tmp_path: Path) -> None:
    result = run(tmp_path, [unique("a", 10)])
    assert not Path(f"{result['prefix']}.cov.plus.bedgraph").exists()


# Edge cases
def test_empty_bam(tmp_path: Path) -> None:
    result = run(tmp_path, [])
    prefix = result["prefix"]
    for name in ("ends.plus", "ends.minus", "cap.plus", "cap.minus"):
        assert Path(f"{prefix}.{name}.bedgraph").read_bytes() == b""
    assert result["stats"]["reads_total"] == 0
    assert not any(result["stats"]["read_fates"].values())


def test_reference_missing_from_fasta_is_an_error(tmp_path: Path) -> None:
    bam, _ = build_inputs(tmp_path, [unique("a", 10)])
    small = tmp_path / "small.fa"
    small.write_text(f">chr1\n{G1}\n")
    pysam.faidx(str(small))
    with pytest.raises(SystemExit, match="missing from FASTA"):
        end_counts.main(
            ["--bam", str(bam), "--fasta", str(small), "--library", "5P", "--prefix", "x"]
        )


def test_fully_soft_clipped_record_is_an_error(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="no aligned bases"):
        run(tmp_path, [sam("r", 0, "chr1", 10, "20S", tags="NH:i:1")])


def test_hard_clip_is_ignored(tmp_path: Path) -> None:
    result = run(tmp_path, [unique("r", 50, cigar="5H20M")])
    assert result[("ends", "+")] == {("chr1", 50): "1"}
    assert result[("cap", "+")] == {}


def test_zero_weight_read_counted_but_not_written(tmp_path: Path) -> None:
    result = run(tmp_path, [sam("r", 0, "chr1", 50, "20M", tags="NH:i:1\tXW:f:0")])
    assert fates(result) == {"unique": 1}
    assert result[("ends", "+")] == {}


def test_fates_sum_to_reads_total(tmp_path: Path) -> None:
    records = [
        *tie_records(["t1", "t2"]),
        unique("short", 50, cigar="10M"),
        sam("unmapped", 4, "*", -1, "*", "A" * 20),
    ]
    stats = run(tmp_path, records)["stats"]
    assert sum(stats["read_fates"].values()) == stats["reads_total"] == 7


def test_compat_not_implemented_yet(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        run(tmp_path, [unique("a", 10)], extra=["--nanopare_compat"])


def test_stats_schema(tmp_path: Path) -> None:
    stats = run(tmp_path, [unique("a", 10)])["stats"]
    assert stats["schema_version"] == 1
    assert stats["parameters"] == {
        "library": "5P",
        "min_aligned": 16,
        "multimapper_mode": "rescue",
        "decimals": 4,
        "nanopare_compat": False,
    }
    json.dumps(stats)
