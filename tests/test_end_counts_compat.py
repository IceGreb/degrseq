"""Unit tests for specs/end-counts.md, `--nanopare_compat`. Each test names the rule or quirk."""

from pathlib import Path

import end_counts_helpers
import pytest
from end_counts_helpers import fates, run, sam
from test_end_counts import G1, first_pos, multi, tie_records, unique

COMPAT = ["--nanopare_compat"]


def compat(tmp_path: Path, records: list[str], library: str = "5P", extra: list[str] | None = None):
    return run(tmp_path, records, library, COMPAT + (extra or []))


def test_compat_requires_rescue_mode(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        compat(tmp_path, [unique("a", 10)], extra=["--multimapper_mode", "uniform"])


def test_compat_rejects_decimals(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        compat(tmp_path, [unique("a", 10)], extra=["--decimals", "4"])


def test_compat_coverage_warns_about_memory(tmp_path: Path, capsys) -> None:
    compat(tmp_path, [unique("a", 10)], extra=["--write_coverage"])
    assert "warning: --write_coverage with --nanopare_compat" in capsys.readouterr().err


# Rule 4, #24
def test_loci_keyed_by_start_only(tmp_path: Path) -> None:
    records = [
        sam("r1", 0, "chr1", 60, "20M", tags="NH:i:2"),
        sam("r1", 0, "chr2", 60, "20M", tags="NH:i:2"),
    ]
    result = compat(tmp_path, records)
    assert fates(result) == {"unique": 1}
    assert result[("ends", "+")] == {("chr1", 60): "1.0"}
    assert result["stats"]["loci"]["loci_key_collisions"] == 1


def test_colliding_spliced_records_merge_positions(tmp_path: Path) -> None:
    # Same start, different introns: one locus, 5' end at the furthest 3' extent on minus.
    records = [
        sam("r1", 16, "chr1", 50, "10M20N10M", tags="NH:i:2"),
        sam("r1", 16, "chr1", 50, "10M40N10M", tags="NH:i:2"),
    ]
    result = compat(tmp_path, records)
    assert result[("ends", "-")] == {("chr1", 109): "1.0"}


# Rule 9
def test_three_prime_only_clip_discarded(tmp_path: Path) -> None:
    seq = G1[60:80] + "GG"
    result = compat(tmp_path, [sam("r1", 0, "chr1", 60, "20M2S", seq, "NH:i:1")])
    assert fates(result) == {"discarded_3prime_clip_only": 1}


# #34
def test_both_end_clip_with_non_g_three_prime_base_discarded(tmp_path: Path) -> None:
    seq = "GG" + G1[60:80] + "AA"
    result = compat(tmp_path, [sam("r1", 0, "chr1", 60, "2S20M2S", seq, "NH:i:1")])
    assert fates(result) == {"discarded_untemplated_non_g": 1}


def test_three_prime_gs_go_to_cap_track(tmp_path: Path) -> None:
    seq = "GG" + G1[60:80] + "GG"
    result = compat(tmp_path, [sam("r1", 0, "chr1", 60, "2S20M2S", seq, "NH:i:1")])
    assert fates(result) == {"unique": 1}
    assert result[("cap", "+")] == {("chr1", p): "1" for p in (58, 59, 80, 81)}


# Rule 22, #8
def test_cap_counted_per_g_at_each_position(tmp_path: Path) -> None:
    seq = "GGG" + G1[50:70]
    result = compat(tmp_path, [sam("r1", 0, "chr1", 50, "3S20M", seq, "NH:i:1")])
    assert result[("cap", "+")] == {("chr1", p): "1" for p in (47, 48, 49)}
    assert result[("ends", "+")] == {("chr1", 50): "1.0"}


# Rule 11, #7
def test_min_aligned_strict_for_unique_reads(tmp_path: Path) -> None:
    result = compat(tmp_path, [unique("len16", 50, cigar="16M"), unique("len17", 100, cigar="17M")])
    assert fates(result) == {"unique": 1, "discarded_too_short": 1}


def test_min_aligned_inclusive_for_multimapper_loci(tmp_path: Path) -> None:
    result = compat(tmp_path, multi("r1", [("chr1", 60, 0, "16M"), ("chr1", 160, 0, "16M")]))
    assert fates(result) == {"multimapped_naive": 1}


# #33
def test_too_short_unique_read_still_adds_cap(tmp_path: Path) -> None:
    seq = "GG" + G1[60:75]
    result = compat(tmp_path, [sam("r1", 0, "chr1", 60, "2S15M", seq, "NH:i:1")])
    assert fates(result) == {"discarded_too_short": 1}
    assert result[("ends", "+")] == {}
    assert result[("cap", "+")] == {("chr1", 58): "1", ("chr1", 59): "1"}


# Rule 20 compat (G9)
def test_tier_two_read_with_one_long_locus_stays_multimapper(tmp_path: Path) -> None:
    result = compat(tmp_path, multi("r1", [("chr1", 60, 0, "20M"), ("chr1", 160, 0, "10M")]))
    assert fates(result) == {"multimapped_naive": 1}
    assert result[("ends", "+")] == {("chr1", 60): "1.0"}


def test_block_takes_last_reads_loci(tmp_path: Path) -> None:
    records = multi("b1", [("chr1", 50), ("chr1", 150)], seq="A" * 20)
    records += multi("b2", [("chr1", 50), ("chr1", 250)], seq="A" * 20)
    result = compat(tmp_path, records)
    assert result[("ends", "+")] == {("chr1", 50): "1.0", ("chr1", 250): "1.0"}


# #36
def test_body_blocks_keyed_by_cigar(tmp_path: Path) -> None:
    records = multi("x", [("chr1", 50), ("chr1", 150)], seq="A" * 20)
    records += multi("y", [("chr1", 50), ("chr2", 150)], seq="C" * 20)
    result = compat(tmp_path, records, library="BODY")
    ends = result[("ends", "+")]
    assert ("chr1", 150) not in ends
    assert ends[("chr2", 150)] == "1.0"


def test_tie_order_follows_bam_order(tmp_path: Path) -> None:
    a = compat(tmp_path / "a", tie_records(["t1", "t2"]))
    b = compat(tmp_path / "b", tie_records(["t2", "t1"]))
    assert a[("ends", "+")] != b[("ends", "+")]


def test_final_pass_uses_coverage_from_earlier_blocks(tmp_path: Path) -> None:
    # No coverage anywhere: a is split equally, then b follows a's coverage at chr1:205.
    records = multi("a", [("chr1", 10), ("chr1", 200)], seq="A" * 20)
    records += multi("b", [("chr1", 205), ("chr2", 100)], seq="C" * 20)
    result = compat(tmp_path, records)
    assert result[("ends", "+")][("chr1", 205)] == "1.0"
    assert ("chr2", 100) not in result[("ends", "+")]


# G2
def test_read_lost_when_only_covered_locus_is_non_g(tmp_path: Path) -> None:
    p = first_pos(lambda q: G1[q - 1] != "A")
    seq = "GA" + G1[p : p + 20]
    records = [unique("u", p)]
    records += multi("m", [("chr1", p, 0, "2S20M"), ("chr1", 200, 0, "22M")], seq=seq)
    result = compat(tmp_path, records)
    assert fates(result) == {"unique": 1, "multimapped_lost_compat": 1}


# #38
def test_cap_zero_rows_written_as_text(tmp_path: Path) -> None:
    records = [unique(f"u{i}", 60) for i in range(200)] + [unique("v", 200)]
    seq = "G" + G1[60:80]
    records += multi("m", [("chr1", 60, 0, "1S20M"), ("chr1", 200, 0, "1S20M")], seq=seq)
    result = compat(tmp_path, records)
    assert result[("cap", "+")] == {("chr1", 59): "1.0", ("chr1", 199): "0.0"}


# #37
def test_only_value_at_position_zero_not_written(tmp_path: Path) -> None:
    result = compat(tmp_path, [unique("a", 0, chrom="chr2"), unique("b", 50)])
    assert result[("ends", "+")] == {("chr1", 50): "1.0"}


# Rule 23 compat
def test_two_decimals_and_merged_rows(tmp_path: Path) -> None:
    records = [
        unique("a", 50),
        unique("b", 51),
        *multi("m", [("chr1", 100), ("chr1", 160), ("chr2", 60)]),
    ]
    result = compat(tmp_path, records)
    text = Path(f"{result['prefix']}.ends.plus.bedgraph").read_text()
    assert text.splitlines()[0] == "chr1\t50\t52\t1.0"
    assert result[("ends", "+")][("chr1", 100)] == "0.33"
    assert result["stats"]["parameters"]["decimals"] == 2


# #39
def test_body_coverage_sums_rounded_strand_values(tmp_path: Path) -> None:
    # 1/3 per strand at chr1:10-28 (the minus records align 19 nt); distinct starts and CIGARs
    # keep loci and blocks apart.
    loci = [("chr1", 10), ("chr1", 100), ("chr2", 150)]
    records = multi("p", [(c, p, 0, "20M") for c, p in loci], seq="A" * 20)
    records += multi("m", [(c, p, 16, "1S19M") for c, p in loci], seq="C" * 20)
    result = compat(tmp_path, records, library="BODY", extra=["--write_coverage"])
    cov = Path(f"{result['prefix']}.cov.bedgraph").read_text().splitlines()
    assert cov[:2] == ["chr1\t10\t29\t0.66", "chr1\t29\t30\t0.33"]


# DL-046
def test_paired_end_rejected(tmp_path: Path) -> None:
    records = [
        sam("r1", 67, "chr1", 50, "20M", tags="NH:i:1"),
        sam("r1", 131, "chr1", 120, "20M", tags="NH:i:1"),
    ]
    with pytest.raises(SystemExit, match="single-end"):
        compat(tmp_path, records)


# Rule 20 compat: multimapper length check counts CIGAR positions, not reclaimed bases
def test_multimapper_length_ignores_reclaimed_bases(tmp_path: Path) -> None:
    seq = G1[57:74]  # 3 clipped bases that match the genome, then 14 aligned
    records = [
        sam("r1", 0, "chr1", 60, "3S14M", seq, "NH:i:2"),
        sam("r1", 256, "chr1", 150, "17M", seq, "NH:i:2"),
    ]
    result = compat(tmp_path, records)
    assert result[("ends", "+")] == {("chr1", 150): "1.0"}


# #24, DL-048: a collision across chromosomes can put coverage past the first reference's end
def test_coverage_past_reference_end_skipped(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(end_counts_helpers.GENOME, "chr3", "A" * 400)
    records = [
        sam("r1", 0, "chr1", 280, "15M", tags="NH:i:2"),
        sam("r1", 0, "chr3", 280, "40M", tags="NH:i:2"),
    ]
    result = compat(tmp_path, records, extra=["--write_coverage"])
    text = Path(f"{result['prefix']}.cov.plus.bedgraph").read_text()
    assert text == "chr1\t280\t300\t1.0\n"
