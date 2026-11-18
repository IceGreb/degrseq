"""Unit tests for specs/kde-peaks.md. Each test names the rule or edge case."""

import json
import math
import random
from pathlib import Path

import numpy as np
import pytest

from degrseq import kde_peaks
from degrseq.kde_peaks import density_peaks, kernel

LENGTHS = "c\t1000\n"
NARROW = ["--bandwidth", "1", "--kernel_extent", "1"]  # W = 1: weights e^-1, 1, e^-1


def run(
    tmp_path: Path,
    ends_plus: str = "",
    ends_minus: str = "",
    body_plus: str = "",
    body_minus: str = "",
    scale: str = "1",
    lengths: str = LENGTHS,
    extra: list[str] | None = None,
) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    files = {
        "ends_plus": ends_plus,
        "ends_minus": ends_minus,
        "body_plus": body_plus,
        "body_minus": body_minus,
        "scale_factor": scale + "\n",
        "lengths": lengths,
    }
    args = []
    for name, text in files.items():
        path = tmp_path / f"{name}.txt"
        path.write_text(text)
        args += [f"--{name}", str(path)]
    prefix = tmp_path / "out"
    kde_peaks.main(args + ["--prefix", str(prefix), "--write_density"] + (extra or []))
    bed = [r.split("\t") for r in Path(f"{prefix}.end_features.bed").read_text().splitlines()]
    dens = {
        label: [
            r.split("\t")
            for r in Path(f"{prefix}.density_{label}.bedgraph").read_text().splitlines()
        ]
        for label in ("plus", "minus")
    }
    stats = json.loads(Path(f"{prefix}.kde-peaks.stats.json").read_text())
    return {"bed": bed, "density": dens, "stats": stats}


def per_nt(rows: list[list[str]]) -> dict[int, float]:
    return {p: float(r[3]) for r in rows for p in range(int(r[1]), int(r[2]))}


# Rule 1: reverse pairs 5P plus with BODY minus; forward with BODY plus; unstranded with both
@pytest.mark.parametrize(
    ("strand", "body", "features"),
    [
        ("reverse", "body_minus", 0),
        ("reverse", "body_plus", 1),
        ("forward", "body_plus", 0),
        ("forward", "body_minus", 1),
        ("unstranded", "body_plus", 0),
        ("unstranded", "body_minus", 0),
    ],
)
def test_pairing(tmp_path: Path, strand: str, body: str, features: int) -> None:
    result = run(
        tmp_path,
        ends_plus="c\t500\t501\t10\n",
        **{body: "c\t500\t501\t10\n"},
        extra=["--body_strand", strand, "--min_rpm", "0"],
    )
    assert len(result["bed"]) == features


# Rules 2, 5: compat rounds the difference to 2 decimals before smoothing; default is exact
def test_difference_rounding(tmp_path: Path) -> None:
    kw = {"ends_plus": "c\t500\t501\t1\n", "scale": "0.3924"}
    exact = per_nt(run(tmp_path / "d", **kw)["density"]["plus"])
    compat = per_nt(run(tmp_path / "c", **kw, extra=["--nanopare_compat"])["density"]["plus"])
    assert exact[500] == round(0.3924 * kernel(15, 3, False)[45], 3)
    assert compat[500] == round(0.39 * kernel(15, 3, True)[45], 3)


# Rule 2, quirk 37: compat drops a track whose only covered position is 0
def test_position_zero_only(tmp_path: Path) -> None:
    kw = {"ends_plus": "c\t0\t1\t10\n", "extra": [*NARROW, "--min_length", "1"]}
    assert len(run(tmp_path / "d", **kw)["bed"]) == 1
    kw["extra"] = [*kw["extra"], "--nanopare_compat"]
    assert run(tmp_path / "c", **kw)["bed"] == []


# Rule 3: weights sum to 1 and are symmetric; compat drops the +W weight
def test_kernel_shape() -> None:
    k = kernel(15, 3, False)
    assert len(k) == 91 and math.isclose(k.sum(), 1) and np.allclose(k, k[::-1])
    c = kernel(15, 3, True)
    assert len(c) == 90 and math.isclose(c[0], k[0], rel_tol=1e-6)


# Rules 3-4: a unit point spreads symmetrically over -45..+45; compat stops at +44
def test_point_spread(tmp_path: Path) -> None:
    d = per_nt(run(tmp_path / "d", ends_plus="c\t500\t501\t100\n")["density"]["plus"])
    assert min(d) == 455 and max(d) == 545 and d[455] == d[545]
    c = per_nt(
        run(tmp_path / "c", ends_plus="c\t500\t501\t100\n", extra=["--nanopare_compat"])["density"][
            "plus"
        ]
    )
    assert min(c) == 455 and max(c) == 544 and (c[455], c[456], c[544]) == (0.174, 0.186, 0.186)


# Rule 4: density stops at the chromosome end in default, not in compat; nothing before 0
def test_chromosome_edges(tmp_path: Path) -> None:
    kw = {"ends_plus": "c\t10\t11\t100\nc\t995\t996\t100\n"}
    d = per_nt(run(tmp_path / "d", **kw)["density"]["plus"])
    assert min(d) == 0 and max(d) == 999
    c = per_nt(run(tmp_path / "c", **kw, extra=["--nanopare_compat"])["density"]["plus"])
    assert min(c) == 0 and max(c) == 1039


# Rule 4: negative differences lower the density
def test_body_lowers_density(tmp_path: Path) -> None:
    alone = per_nt(run(tmp_path / "a", ends_plus="c\t500\t501\t10\n")["density"]["plus"])
    sub = per_nt(
        run(tmp_path / "b", ends_plus="c\t500\t501\t10\n", body_minus="c\t510\t511\t10\n")[
            "density"
        ]["plus"]
    )
    assert sub[500] < alone[500] and max(sub) < max(alone)


# Rule 6: length >= min_length kept, one shorter dropped; default end exclusive, compat inclusive
def test_feature_length_and_end(tmp_path: Path) -> None:
    eight = "".join(f"c\t{p}\t{p + 1}\t1\n" for p in range(100, 108))  # positive 99..108 = 10 nt
    seven = "".join(f"c\t{p}\t{p + 1}\t1\n" for p in range(300, 307))  # positive 299..307 = 9 nt
    d = run(tmp_path / "d", ends_plus=eight + seven, extra=NARROW)["bed"]
    assert [(r[1], r[2]) for r in d] == [("99", "109")]
    # compat drops the +W weight (quirk 10): the same reads give 99..107 (9 nt) and 299..306 (8 nt)
    c = run(
        tmp_path / "c",
        ends_plus=eight + seven,
        extra=[*NARROW, "--min_length", "9", "--nanopare_compat"],
    )["bed"]
    assert [(r[1], r[2]) for r in c] == [("99", "107")]


# Rule 7: run entered by a rise and followed by a fall; ties by larger position
@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([1, 1, 3, 1, 1, 3, 1, 1, 1, 1], [5, 2]),
        ([2, 2, 2, 5, 5, 5, 1, 1, 1, 1], [5]),
        ([1, 3, 3, 1, 3, 3], [2]),
        ([3, 2, 1], [0]),
        ([1, 1, 1], [2]),
        ([1, 2, 3], [2]),
    ],
)
def test_density_peaks(values: list[int], expected: list[int]) -> None:
    assert density_peaks([float(v) for v in values]) == expected


# Rule 7: column 9 holds at most 9 further peaks
def test_other_peaks_capped(tmp_path: Path) -> None:
    ends = "".join(f"c\t{p}\t{p + 1}\t5\n" for p in range(100, 160, 2))  # 30 equal peaks
    row = run(tmp_path, ends_plus=ends, extra=[*NARROW, "--min_length", "1"])["bed"][0]
    assert row[7] == "59" and row[8] == "57,55,53,51,49,47,45,43,41"


# Rule 8: dominant 5P position, ties to the larger position
def test_dominant_tie(tmp_path: Path) -> None:
    ends = "c\t100\t101\t4\nc\t102\t103\t4\n"
    row = run(tmp_path, ends_plus=ends, extra=[*NARROW, "--min_length", "1"])["bed"][0]
    assert (row[1], row[4]) == ("99", "3")


# Rules 8-10, quirk 12: compat ignores the feature's last nt, so a read there gives RPM 0
def test_read_on_last_nt(tmp_path: Path) -> None:
    kw = {"ends_plus": "c\t100\t101\t10\n", "body_minus": "c\t102\t103\t50\n"}
    d = run(tmp_path / "d", **kw, extra=[*NARROW, "--min_length", "1"])["bed"]
    assert [(r[1], r[2], r[4]) for r in d] == [("99", "101", "1")]
    c = run(tmp_path / "c", **kw, extra=[*NARROW, "--min_length", "1", "--nanopare_compat"])
    assert c["bed"] == [] and c["stats"]["plus"]["features_without_reads"] == 1


# Rule 9: default N = plus + minus; compat N = int(2 x plus)
def test_rpm_denominator(tmp_path: Path) -> None:
    kw = {"ends_plus": "c\t100\t101\t3\n", "ends_minus": "c\t500\t501\t1\n", "extra": []}
    d = run(tmp_path / "d", **kw)
    assert d["stats"]["total_5p"] == 4 and d["bed"][0][6] == "750000.00"
    c = run(tmp_path / "c", **{**kw, "extra": [*kw["extra"], "--nanopare_compat"]})
    assert c["stats"]["total_5p"] == 6 and c["bed"][0][6] == "500000.0"


# Rules 6, 10: names count features before the RPM filter
def test_names_count_filtered(tmp_path: Path) -> None:
    ends = "c\t100\t101\t1\nc\t500\t501\t1000\n"
    bed = run(tmp_path, ends_plus=ends, extra=[*NARROW, "--min_length", "1", "--min_rpm", "10000"])[
        "bed"
    ]
    assert [r[3] for r in bed] == ["5P.plus.c.2"]


# Rule 11: chromosomes in --lengths order; equal starts plus first
def test_order(tmp_path: Path) -> None:
    lengths = "z\t1000\na\t1000\n"
    one = "z\t100\t101\t1\na\t100\t101\t1\n"
    bed = run(
        tmp_path,
        ends_plus=one,
        ends_minus=one,
        lengths=lengths,
        extra=[*NARROW, "--min_length", "1"],
    )["bed"]
    assert [(r[0], r[5]) for r in bed] == [("z", "+"), ("z", "-"), ("a", "+"), ("a", "-")]


# Rule 12: default fixed decimals; compat repr
def test_number_format(tmp_path: Path) -> None:
    d = run(tmp_path / "d", ends_plus="c\t500\t501\t1\n")["density"]["plus"]
    assert all(len(r[3].split(".")[1]) == 3 for r in d)
    c = run(tmp_path / "c", ends_plus="c\t500\t501\t1\n", extra=["--nanopare_compat"])["density"][
        "plus"
    ]
    assert any(len(r[3].split(".")[1]) < 3 for r in c)


# Rule 13: row order does not matter
def test_row_order(tmp_path: Path) -> None:
    rows = [f"c\t{p}\t{p + 1}\t{1 + p % 5}\n" for p in range(100, 400, 7)]
    a = run(tmp_path / "a", ends_plus="".join(rows))
    random.Random(1).shuffle(rows)
    b = run(tmp_path / "b", ends_plus="".join(rows))
    assert a == b


# Edge cases: invalid input is an error
@pytest.mark.parametrize(
    "kw",
    [
        {"ends_plus": "x\t1\t2\t1\n"},
        {"ends_plus": "c\t999\t1001\t1\n"},
        {"body_plus": "c\t1\t2\t-1\n"},
        {"ends_plus": "c\t1\t2\tnan\n"},
        {"scale": "0"},
    ],
)
def test_invalid_input(tmp_path: Path, kw: dict) -> None:
    with pytest.raises(SystemExit):
        run(tmp_path, **kw)


# Edge case: no 5P values -> empty output, not an error; compat divides by zero -> error
def test_no_5p(tmp_path: Path) -> None:
    assert run(tmp_path / "a", body_plus="c\t1\t2\t3\n")["bed"] == []
    with pytest.raises(SystemExit):
        run(tmp_path / "b", ends_minus="c\t500\t501\t1\n", extra=["--nanopare_compat"])


# Edge case: chromosome names are literal strings (DL-064)
def test_regex_chromosome_name(tmp_path: Path) -> None:
    bed = run(tmp_path, ends_plus="Su(var)3-9\t500\t501\t5\n", lengths="Su(var)3-9\t1000\n")["bed"]
    assert bed[0][0] == "Su(var)3-9" and bed[0][3] == "5P.plus.Su(var)3-9.1"


# Parameters: bandwidth above the original's cap of 30 is an error
def test_bandwidth_cap(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        run(tmp_path, ends_plus="c\t500\t501\t5\n", extra=["--bandwidth", "31"])
