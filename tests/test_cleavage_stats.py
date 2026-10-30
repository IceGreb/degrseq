"""Unit tests for specs/cleavage-stats.md. Each test names the rule or edge case."""

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from degrseq import cleavage_stats
from degrseq.cleavage_stats import Site

HEADER = (
    "sRNA\tTranscript\tSlice.site\tAllenScore\tflankMax.rpm\tslice.site.rpm\tslice.site.hit.norm\n"
)
SRNAS = "sRNA\tclass\tfamily\nmiR156a\t\t\nmiR156b\t\t\nmiR2\t\t\nTAS3A\t\t\n"


def table(rows: list[tuple]) -> str:
    return HEADER + "".join("\t".join(map(str, r)) + "\t1\n" for r in rows)


def null(srna: str, n: int = 4, fc_rpm: float = 0, allen: float = 6) -> list[tuple]:
    """n shuffled rows for srna, each with fold change (fc_rpm + 1) / 1."""
    return [(srna, "TX9.1", 10 * i, allen, 0, fc_rpm) for i in range(n)]


def run(
    tmp_path: Path,
    real: list[tuple],
    shuffled: list[tuple],
    real_w50: list[tuple] | None = None,
    srnas: str = SRNAS,
    extra: list[str] | None = None,
    files: dict[str, str] | None = None,
) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    paths = {
        "real_sites_w20": table(real),
        "real_sites_w50": table(real if real_w50 is None else real_w50),
        "shuffled_sites_w20": table(shuffled),
        "shuffled_sites_w50": table(shuffled),
        "srna_table": srnas,
        **(files or {}),
    }
    argv = []
    for flag, text in paths.items():
        path = tmp_path / f"{flag}.tsv"
        path.write_text(text)
        argv += [f"--{flag}", str(path)]
    prefix = tmp_path / "out"
    cleavage_stats.main([*argv, "--prefix", str(prefix), *(extra or [])])
    out = {"stats": json.loads(Path(f"{prefix}.cleavage-stats.stats.json").read_text())}
    for cls in cleavage_stats.CLASSES:
        for kind in ("detected", "significant"):
            text = Path(f"{prefix}.{cls}.{kind}_sites.tsv").read_text()
            lines = text.splitlines()[1:]
            out[f"{cls}.{kind}"] = [
                dict(zip(cleavage_stats.HEADER, line.split("\t"), strict=True)) for line in lines
            ]
    return out


# Rule 1
def test_fold_change_pseudocount(tmp_path: Path) -> None:
    out = run(tmp_path, [("miR156a", "AT1G01.1", 100, 2, 1, 3)], null("miR156a"))
    assert out["miRNA.detected"][0]["fold.change"] == "2"


# Rules 2, 13
def test_class_from_table_or_name(tmp_path: Path) -> None:
    srnas = "sRNA\tclass\tfamily\nmiR156a\t\t\nTAS3A\t\t\name-miR-1\t\t\nlet-7\tmiRNA\tlet7\n"
    real = [
        (s, "AT1G01.1", 100 * i, 2, 0, 5)
        for i, s in enumerate(["miR156a", "TAS3A", "ame-miR-1", "let-7", "other"])
    ]
    shuffled = [r for s in ("miR156a", "TAS3A", "let-7") for r in null(s)]
    out = run(tmp_path, real, shuffled, srnas=srnas)
    assert sorted(r["sRNA"] for r in out["miRNA.detected"]) == ["let-7", "miR156a"]
    assert [r["sRNA"] for r in out["tasiRNA.detected"]] == ["TAS3A"]
    assert out["stats"]["n_srnas_unclassified"] == 1


@pytest.mark.parametrize(
    ("srnas", "message"),
    [
        ("sRNA\nmiR1\nmiR1\n", "duplicate sRNA miR1"),
        ("sRNA\tclass\nmiR1\tsnoRNA\n", "unknown class 'snoRNA'"),
        ("name\nmiR1\n", "missing column"),
    ],
)
def test_bad_srna_table_exits(tmp_path: Path, srnas: str, message: str) -> None:
    with pytest.raises(SystemExit, match=message):
        run(tmp_path, [], [], srnas=srnas)


# Rule 3
def test_self_self_removed_and_missing_precursor_kept(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 0, 5), ("miR2", "AT1G01.2", 100, 2, 0, 5)]
    files = {"precursors": "sRNA\tprecursor\nmiR156a\tAT1G01\n"}
    out = run(tmp_path, real, null("miR156a") + null("miR2"), files=files)
    assert [r["sRNA"] for r in out["miRNA.detected"]] == ["miR2"]
    assert out["stats"]["classes"]["miRNA"]["windows"]["20"]["n_after_self_self"] == 1


def test_several_precursors_per_srna(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 0, 5), ("miR156a", "AT1G02.1", 100, 2, 0, 5)]
    files = {"precursors": "sRNA\tprecursor\nmiR156a\tAT1G01\nmiR156a\tAT1G02\n"}
    out = run(tmp_path, real, null("miR156a"), files=files)
    assert out["miRNA.detected"] == []


def test_srnas_not_in_table_warned(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = run(tmp_path, [("ath-miR156a", "AT1G01.1", 100, 2, 0, 5)], null("ath-miR156a"))
    assert out["stats"]["n_srnas_not_in_table"] == 1
    assert "1 sRNA(s) with real sites are not in --srna_table" in capsys.readouterr().err


def test_negative_collapse_nt_rejected(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        run(tmp_path, [], [], extra=["--collapse_nt", "-1"])


# Rule 4
@pytest.mark.parametrize(
    ("rows", "kept"),
    [
        ([("miR2", 103, 3, 9), ("miR156a", 100, 2, 1)], {"miR156a"}),  # lower Allen wins
        ([("miR2", 103, 2, 9), ("miR156a", 100, 2, 1)], {"miR2"}),  # then higher rpm
        ([("miR2", 103, 2, 1), ("miR156a", 100, 2, 1)], {"miR156a"}),  # then sRNA name
        ([("miR2", 105, 3, 1), ("miR156a", 100, 2, 1)], {"miR156a"}),  # 5 nt apart: redundant
        ([("miR2", 106, 3, 1), ("miR156a", 100, 2, 1)], {"miR156a", "miR2"}),  # 6 nt: kept
    ],
)
def test_collapse_order(tmp_path: Path, rows: list[tuple], kept: set[str]) -> None:
    real = [(s, "AT1G01.1", site, allen, 0, rpm) for s, site, allen, rpm in rows]
    out = run(tmp_path, real, null("miR156a") + null("miR2"))
    assert {r["sRNA"] for r in out["miRNA.detected"]} == kept


def test_collapse_is_per_transcript() -> None:
    sites = [Site("a", "T1", 100, 2, 1, 1), Site("b", "T2", 100, 2, 1, 1)]
    assert len(cleavage_stats.collapse(sites, 5)) == 2


# Rule 6
def test_ecdf_strict_for_fc_inclusive_for_allen() -> None:
    real = [Site("s", "T", 1, 5, 0, 2.0)]
    null = {"s": (np.array([1.0, 2, 2, 3]), np.array([4.0, 5, 5, 6]))}
    cleavage_stats.window_pvalues(real, null, set())
    # R: 1 - ecdf(c(1,2,2,3))(2) = 0.25; ecdf(c(4,5,5,6))(5) = 0.75; Fisher = 0.50137058129468837
    assert real[0].adj == pytest.approx(0.50137058129468837, rel=1e-14)


# Rule 7
def test_fisher_zero_gives_zero() -> None:
    real = [Site("s", "T", 1, 2, 0, 9.0)]
    cleavage_stats.window_pvalues(real, {"s": (np.array([1.0]), np.array([3.0]))}, set())
    assert real[0].adj == 0


# Rule 8
def test_bh_matches_p_adjust() -> None:
    p = np.array([0.001, 0.2, 0.05, 0.05, 0.9, 0])
    expected = [0.003, 0.24, 0.075, 0.075, 0.9, 0]  # R p.adjust(method = "fdr")
    assert cleavage_stats.bh(p).tolist() == pytest.approx(expected, rel=1e-15)


# Rule 6, D3
def test_empty_null_gives_p_one_and_warns(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = run(tmp_path, [("miR156a", "AT1G01.1", 100, 2, 0, 5)], null("miR2"))
    assert out["miRNA.detected"][0]["adj.p.val"] == "1"
    assert out["stats"]["srnas_empty_null"] == ["miR156a"]
    assert "miR156a has real sites but no shuffled sites" in capsys.readouterr().err


# Rule 9
def test_significant_filter_and_thresholds(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 0, 5), ("miR2", "AT1G02.1", 100, 2, 0, 0.05)]
    shuffled = null("miR156a") + null("miR2")
    out = run(tmp_path, real, shuffled)
    assert [r["sRNA"] for r in out["miRNA.significant"]] == ["miR156a"]
    assert out["stats"]["classes"]["miRNA"]["n_significant"] == 1
    out = run(tmp_path / "strict", real, shuffled, extra=["--fc_min", "10"])
    assert out["miRNA.significant"] == []


# Rule 10
def test_window_join_takes_smaller_adj_p(tmp_path: Path) -> None:
    w20 = [("miR156a", "AT1G01.1", 100, 2, 0, 5)]
    w50 = [("miR156a", "AT1G01.1", 100, 2, 0, 0.5)]
    shuffled = null("miR156a", fc_rpm=1, allen=1)
    out = run(tmp_path, w20, shuffled, real_w50=w50)
    assert out["miRNA.detected"][0]["slice.site.rpm"] == "5"
    out = run(tmp_path / "swap", w50, shuffled, real_w50=w20)
    assert out["miRNA.detected"][0]["slice.site.rpm"] == "5"


def test_window_tie_takes_w20(tmp_path: Path) -> None:
    w20 = [("miR156a", "AT1G01.1", 100, 2, 0, 5)]
    w50 = [("miR156a", "AT1G01.1", 100, 2, 0, 6)]
    out = run(tmp_path, w20, null("miR156a"), real_w50=w50)
    assert out["miRNA.detected"][0]["adj.p.val"] == "0"
    assert out["miRNA.detected"][0]["slice.site.rpm"] == "5"


# Edge case, Q10
def test_site_in_w50_only(tmp_path: Path) -> None:
    out = run(tmp_path, [], null("miR156a"), real_w50=[("miR156a", "AT1G01.1", 100, 2, 0, 5)])
    assert out["miRNA.detected"][0]["target.transcript"] == "AT1G01.1"


def test_class_without_sites_writes_headers(tmp_path: Path) -> None:
    out = run(tmp_path, [("miR156a", "AT1G01.1", 100, 2, 0, 5)], null("miR156a"))
    assert out["tasiRNA.detected"] == [] and out["tasiRNA.significant"] == []


# Rule 11
LABELS = {
    "gene_aliases": "locus_name\tsymbol\nAT1G01\tSPL1\nAT1G01\tSPL1B\n",
    "validated_targets": "sRNA\ttargetAGI\nmiR156\tAT1G01\n",
}


def test_labels(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 0, 5), ("miR2", "AT1G02.1", 100, 2, 0, 5)]
    out = run(tmp_path, real, null("miR156a") + null("miR2"), files=LABELS)
    rows = {r["sRNA"]: r for r in out["miRNA.detected"]}
    assert (rows["miR156a"]["verified."], rows["miR156a"]["common.name"]) == ("yes", "SPL1")
    assert rows["miR156a"]["interaction"] == "miR156a:SPL1"
    assert (rows["miR2"]["verified."], rows["miR2"]["common.name"]) == ("no", "NA")
    assert rows["miR2"]["interaction"] == "miR2:AT1G02.1"
    assert out["stats"]["classes"]["miRNA"]["n_significant_validated"] == 1


def test_validated_family_from_table(tmp_path: Path) -> None:
    srnas = "sRNA\tclass\tfamily\nmiR2934-3p\t\tmiR2934\n"
    files = {"validated_targets": "sRNA\ttargetAGI\nmiR2934\tAT1G01\n"}
    out = run(
        tmp_path,
        [("miR2934-3p", "AT1G01.1", 100, 2, 0, 5)],
        null("miR2934-3p"),
        srnas=srnas,
        files=files,
    )
    assert out["miRNA.detected"][0]["verified."] == "yes"


# Rule 11, D2 / Q3
@pytest.mark.parametrize(("compat", "expected"), [(False, "SPL1"), (True, "NA")])
def test_one_window_site_common_name(tmp_path: Path, compat: bool, expected: str) -> None:
    extra = ["--nanopare_compat"] if compat else []
    out = run(
        tmp_path,
        [("miR156a", "AT1G01.1", 100, 2, 0, 5)],
        null("miR156a"),
        real_w50=[],
        files=LABELS,
        extra=extra,
    )
    assert out["miRNA.detected"][0]["common.name"] == expected


@pytest.mark.parametrize(("compat", "expected"), [(False, "yes"), (True, "no")])
def test_verified_when_windows_disagree(tmp_path: Path, compat: bool, expected: str) -> None:
    w20 = [("miR2", "AT1G01.1", 100, 2, 0, 5)]
    w50 = [("miR156a", "AT1G01.1", 100, 2, 0, 5)]
    extra = ["--nanopare_compat"] if compat else []
    out = run(
        tmp_path, w20, null("miR156a") + null("miR2"), real_w50=w50, files=LABELS, extra=extra
    )
    assert out["miRNA.detected"][0]["verified."] == expected


# Rule 11, Q9
@pytest.mark.parametrize(("compat", "expected"), [(False, "miR156a:SPL1"), (True, "miR2:SPL1")])
def test_interaction_srna(tmp_path: Path, compat: bool, expected: str) -> None:
    w20 = [("miR2", "AT1G01.1", 100, 2, 0, 5)]
    w50 = [("miR156a", "AT1G01.1", 100, 2, 0, 5)]
    shuffled = null("miR156a") + null("miR2", fc_rpm=9, allen=1)
    extra = ["--nanopare_compat"] if compat else []
    out = run(tmp_path, w20, shuffled, real_w50=w50, files=LABELS, extra=extra)
    row = out["miRNA.detected"][0]
    assert row["sRNA"] == "miR156a"
    assert row["interaction"] == expected


# Rule 12
@pytest.mark.parametrize(
    ("compat", "expected"),
    [(False, ["AT1G01.1_9", "AT1G01.1_100"]), (True, ["AT1G01.1_100", "AT1G01.1_9"])],
)
def test_output_order(tmp_path: Path, compat: bool, expected: list[str]) -> None:
    # Both adj p = 0; default breaks the tie on fc, compat on the text Transcript_Slice.site.
    real = [("miR2", "AT1G01.1", 100, 2, 0, 5), ("miR156a", "AT1G01.1", 9, 2, 0, 9)]
    extra = ["--nanopare_compat"] if compat else []
    out = run(tmp_path, real, null("miR156a") + null("miR2"), extra=extra)
    order = [f"{r['target.transcript']}_{r['slice.site']}" for r in out["miRNA.detected"]]
    assert order == expected


# Rule 14
@pytest.mark.parametrize(
    ("value", "text"),
    [
        (0.0001, "1e-04"),
        (1e-5, "1e-05"),
        (0.1 + 0.2, "0.3"),
        (123456789012345.0, "123456789012345"),
        (1 / 3, "0.333333333333333"),
        (2.5, "2.5"),
        (3.0, "3"),
        (0.0, "0"),
        (14.5175, "14.5175"),
        (0.000123456, "0.000123456"),
        (1e15, "1e+15"),
        (100000.0, "1e+05"),
        (1234567.1, "1234567.1"),
    ],
)
def test_r_number_format(value: float, text: str) -> None:
    assert cleavage_stats.r_number(value) == text  # reference: R 4.5.3 write.table


# Rules 1, 5: input handling
def test_missing_shuffled_values_become_zero(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 0, 0.5)]
    shuffled = [("miR156a", "TX", 1, 6, "NA", 5), ("miR156a", "TX", 2, "", 0, 0)]
    out = run(tmp_path, real, shuffled)
    assert out["miRNA.detected"][0]["adj.p.val"] == "0"


def test_nan_and_missing_slice_in_shuffled_rows(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 0, 0.5)]
    shuffled = [("miR156a", "TX", "", "nan", 0, 5), ("miR156a", "TX", "NA", 1, 0, 0)]
    out = run(tmp_path, real, shuffled)
    assert out["stats"]["classes"]["miRNA"]["windows"]["20"]["n_shuffled_rows"] == 2


@pytest.mark.parametrize("value", ["inf", "-1", "NaN"])
def test_bad_real_value_exits(tmp_path: Path, value: str) -> None:
    with pytest.raises(SystemExit, match="flankMax"):
        run(tmp_path, [("miR156a", "AT1G01.1", 100, 2, value, 5)], null("miR156a"))


def test_corrupt_gzip_exits(tmp_path: Path) -> None:
    tmp_path.mkdir(parents=True, exist_ok=True)
    bad = tmp_path / "bad.tsv.gz"
    bad.write_text("not gzip")
    with pytest.raises(SystemExit, match=r"bad\.tsv\.gz"):
        run(tmp_path, [], [], extra=["--shuffled_sites_w20", str(bad)])


def test_missing_real_value_exits(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match=r"missing slice\.site\.rpm in a real-site row"):
        run(tmp_path, [("miR156a", "AT1G01.1", 100, 2, 0, "NA")], null("miR156a"))


def test_gzip_input(tmp_path: Path) -> None:
    path = tmp_path / "null.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(table(null("miR156a")))
    real = tmp_path / "real.tsv"
    real.write_text(table([("miR156a", "AT1G01.1", 100, 2, 0, 5)]))
    srnas = tmp_path / "srnas.tsv"
    srnas.write_text(SRNAS)
    prefix = tmp_path / "out"
    cleavage_stats.main(
        ["--real_sites_w20", str(real), "--real_sites_w50", str(real),
         "--shuffled_sites_w20", str(path), "--shuffled_sites_w50", str(path),
         "--srna_table", str(srnas), "--prefix", str(prefix)]
    )  # fmt: skip
    stats = json.loads(Path(f"{prefix}.cleavage-stats.stats.json").read_text())
    assert stats["classes"]["miRNA"]["windows"]["20"]["n_shuffled_rows"] == 4


def test_output_independent_of_hash_seed(tmp_path: Path) -> None:
    real = [
        (s, f"AT1G0{i}.1", 100, 2, 0, 5)
        for i, s in enumerate(["miR156a", "miR156b", "miR2", "TAS3A"])
    ]
    shuffled = [r for s in ("miR156a", "miR156b", "miR2", "TAS3A") for r in null(s)]
    paths = {}
    for name, text in {"real": table(real), "null": table(shuffled), "srnas": SRNAS}.items():
        paths[name] = tmp_path / f"{name}.tsv"
        paths[name].write_text(text)
    outputs = []
    for hash_seed in ("0", "123"):
        prefix = tmp_path / f"h{hash_seed}"
        subprocess.run(
            [sys.executable, "-m", "degrseq", "cleavage-stats",
             "--real_sites_w20", str(paths["real"]), "--real_sites_w50", str(paths["real"]),
             "--shuffled_sites_w20", str(paths["null"]), "--shuffled_sites_w50", str(paths["null"]),
             "--srna_table", str(paths["srnas"]), "--prefix", str(prefix)],
            check=True, env={**os.environ, "PYTHONHASHSEED": hash_seed},
        )  # fmt: skip
        suffixes = (
            ".miRNA.detected_sites.tsv",
            ".tasiRNA.detected_sites.tsv",
            ".cleavage-stats.stats.json",
        )
        outputs.append([Path(f"{prefix}{suffix}").read_bytes() for suffix in suffixes])
    assert outputs[0] == outputs[1]


# Rule 15
def test_window_tests_and_shuffled_values(tmp_path: Path) -> None:
    real = [("miR156a", "AT1G01.1", 100, 2, 1, 5)]
    shuffled = null("miR156a", n=3) + null("miR2", n=2)
    run(tmp_path, real, shuffled)
    tests = (tmp_path / "out.miRNA.window_tests.tsv").read_text().splitlines()
    assert tests[0].split("\t") == list(cleavage_stats.WINDOW_HEADER)
    assert [r.split("\t")[:2] + r.split("\t")[7:9] for r in tests[1:]] == [
        ["20", "miR156a", "3", "3"],
        ["50", "miR156a", "3", "3"],
    ]
    values = (tmp_path / "out.miRNA.shuffled_values.tsv").read_text().splitlines()
    assert values == ["window\tsRNA\tfold.change\tAllenScore"] + [
        f"{w}\tmiR156a\t1\t6" for w in (20, 20, 20, 50, 50, 50)
    ]
