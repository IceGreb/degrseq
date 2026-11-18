"""Unit tests for specs/cap-class.md: one per rule and edge case, handmade fixtures."""

import json
import random
from pathlib import Path

import pytest

from degrseq import cap_class


def write(path: Path, rows: list[tuple]) -> str:
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def run(
    tmp_path: Path,
    features: list[tuple],
    plus: list[list[tuple]],
    minus: list[list[tuple]] | None = None,
    ug_plus: list[list[tuple]] | None = None,
    ug_minus: list[list[tuple]] | None = None,
    extra: list[str] | None = None,
    lengths: tuple = (("c1", 1000),),
    name: str = "out",
) -> dict:
    n = len(plus)
    minus = minus or [[] for _ in range(n)]
    ug_plus = ug_plus or [[] for _ in range(n)]
    ug_minus = ug_minus or [[] for _ in range(n)]
    args = ["--lengths", write(tmp_path / "len.tsv", list(lengths))]
    args += ["--features", write(tmp_path / "f.bed", features)]
    for flag, tracks in (
        ("--ends_plus", plus),
        ("--ends_minus", minus),
        ("--ug_plus", ug_plus),
        ("--ug_minus", ug_minus),
    ):
        args += [flag] + [write(tmp_path / f"{flag[2:]}{i}.bg", t) for i, t in enumerate(tracks)]
    prefix = tmp_path / name
    cap_class.main(args + ["--prefix", str(prefix)] + (extra or []))
    out = {}
    for f in (
        "plus.bedgraph",
        "minus.bedgraph",
        "uG.plus.bedgraph",
        "peaks.bed",
        "capped.bed",
        "noncapped.bed",
    ):
        out[f] = [line.split("\t") for line in Path(f"{prefix}.{f}").read_text().splitlines()]
    out["stats"] = json.loads(Path(f"{prefix}.cap-class.stats.json").read_text())
    return out


REGION = [("c1", 100, 120, "g.1", 2, "+")]


def test_rule1_2_sums_default_exact_and_formatted(tmp_path: Path) -> None:
    reps = [
        [("c1", 10, 12, 0.333)],
        [("c1", 10, 11, 0.334), ("c1", 10, 11, 1)],
        [("c1", 11, 12, 0.334)],
    ]
    out = run(tmp_path, REGION, reps)
    assert out["plus.bedgraph"] == [["c1", "10", "11", "1.6670"], ["c1", "11", "12", "0.6670"]]


def test_rule1_2_sums_compat_rounded_repr(tmp_path: Path) -> None:
    reps = [[("c1", 10, 12, 0.333)], [("c1", 10, 12, 0.334)]]
    out = run(tmp_path, REGION, reps, extra=["--nanopare_compat"])
    assert out["plus.bedgraph"] == [["c1", "10", "12", "0.67"]]


def test_rule2_order_default_lengths_compat_sorted(tmp_path: Path) -> None:
    lengths = (("zz", 100), ("aa", 100))
    reps = [[("zz", 5, 6, 1), ("aa", 5, 6, 1)]]
    feats = [("zz", 0, 10, "g.1", 2, "+")]
    d = run(tmp_path, feats, reps, lengths=lengths, name="d")
    c = run(tmp_path, feats, reps, lengths=lengths, name="c", extra=["--nanopare_compat"])
    assert [r[0] for r in d["plus.bedgraph"]] == ["zz", "aa"]
    assert [r[0] for r in c["plus.bedgraph"]] == ["aa", "zz"]


def test_rule3_4_trim_and_rightmost_peak(tmp_path: Path) -> None:
    reps = [[("c1", 105, 106, 3), ("c1", 110, 111, 5), ("c1", 112, 113, 5)]]
    out = run(tmp_path, REGION, reps)
    assert out["peaks.bed"] == [["c1", "105", "113", "g.1", "7", "+"]]
    assert out["capped.bed"] == [] and out["noncapped.bed"][0][1:3] == ["105", "113"]


def test_rule3_zero_signal_region_dropped(tmp_path: Path) -> None:
    out = run(tmp_path, [*REGION, ("c1", 300, 320, "g.2", 2, "+")], [[("c1", 105, 106, 3)]])
    assert len(out["peaks.bed"]) == 1 and out["stats"]["dropped_zero_signal"] == 1


def test_rule4_peak_rightmost_on_minus_strand_too(tmp_path: Path) -> None:
    feats = [("c1", 100, 120, "g.1", 2, "-")]
    out = run(tmp_path, feats, [[]], minus=[[("c1", 101, 102, 4), ("c1", 108, 109, 4)]])
    assert out["peaks.bed"][0][1:6] == ["101", "109", "g.1", "7", "-"]


def test_rule5_plus_strand_cap_at_start_default_vs_compat(tmp_path: Path) -> None:
    five = [[("c1", 105, 106, 10)]]
    d = run(
        tmp_path, REGION, five, ug_plus=[[("c1", 105, 106, 10)]], name="d"
    )  # default track: at 5' end
    c = run(
        tmp_path,
        REGION,
        five,
        ug_plus=[[("c1", 104, 105, 10)]],
        name="c",
        extra=["--nanopare_compat"],
    )
    assert d["capped.bed"][0][1:3] == ["105", "106"] and d["capped.bed"][0][9] == "1.000"
    assert (
        c["noncapped.bed"][0][1:3] == ["105", "107"] and c["noncapped.bed"][0][9] == "0.0"
    )  # end + 1


def test_rule5_compat_minus_strand_g_inside_buffer(tmp_path: Path) -> None:
    feats = [("c1", 100, 120, "g.1", 2, "-")]
    c = run(
        tmp_path,
        feats,
        [[]],
        minus=[[("c1", 110, 111, 4)]],
        ug_minus=[[("c1", 111, 112, 4)]],
        extra=["--nanopare_compat"],
    )
    assert c["capped.bed"][0][1:6] == ["110", "112", "5P.capped.1", "1", "-"]  # end + 1, peak + 1


def test_rule8_cutoff_inclusive(tmp_path: Path) -> None:
    out = run(tmp_path, REGION, [[("c1", 105, 106, 10)]], ug_plus=[[("c1", 105, 106, 1)]])
    assert out["capped.bed"][0][9] == "0.100"


def test_rule8_min_uug_zero_default_all_capped(tmp_path: Path) -> None:
    out = run(tmp_path, REGION, [[("c1", 105, 106, 10)]], extra=["--min_uug", "0"])
    assert len(out["capped.bed"]) == 1 and out["stats"]["min_uug"] == 0


def test_rule8_min_uug_zero_compat_automatic_tie_mean(tmp_path: Path) -> None:
    feats, five, ug, i = [], [], [], 0
    for b in range(25):
        for _ in range(1 if b in (3, 7) else 2):
            s = i * 10
            i += 1
            feats.append(("c1", s, s + 1, f"g.{i}", 2, "+"))
            five.append(("c1", s, s + 1, 100))
            if b:
                ug.append(("c1", s, s + 1, b))
    out = run(tmp_path, feats, [five], ug_plus=[ug], extra=["--nanopare_compat", "--min_uug", "0"])
    assert out["stats"]["min_uug"] == pytest.approx(
        0.05
    )  # bins 3 and 7 tie (same as the original's helper)
    assert len(out["capped.bed"]) == 39


def test_rule7_below_min_reads_dropped(tmp_path: Path) -> None:
    reps = [[("c1", 105, 106, 0.3)], [("c1", 105, 106, 0.3)], [("c1", 105, 106, 0.3)]]
    out = run(tmp_path, REGION, reps)
    assert (
        out["capped.bed"] == out["noncapped.bed"] == []
        and out["stats"]["dropped_below_min_reads"] == 1
    )
    assert len(out["peaks.bed"]) == 1


def test_rule9_same_start_both_strands_are_two_features(tmp_path: Path) -> None:
    feats = [("c1", 100, 120, "g.1", 2, "+"), ("c1", 100, 120, "g.2", 2, "-")]
    out = run(tmp_path, feats, [[("c1", 100, 101, 2)]], minus=[[("c1", 100, 101, 2)]])
    assert [r[5] for r in out["noncapped.bed"]] == ["+", "-"]
    assert [r[3] for r in out["noncapped.bed"]] == ["5P.noncapped.1", "5P.noncapped.2"]


def test_quirk37_position_zero_only(tmp_path: Path) -> None:
    feats = [("c1", 0, 10, "g.1", 2, "+")]
    d = run(tmp_path, feats, [[("c1", 0, 1, 5)]], name="d")
    c = run(tmp_path, feats, [[("c1", 0, 1, 5)]], name="c", extra=["--nanopare_compat"])
    assert d["plus.bedgraph"] == [["c1", "0", "1", "5.0000"]] and len(d["noncapped.bed"]) == 1
    assert c["plus.bedgraph"] == [] and c["stats"]["dropped_zero_signal"] == 1


def test_ug_exceeding_5p_errors_in_default_only(tmp_path: Path, capsys) -> None:
    five, ug = [[("c1", 105, 106, 1)]], [[("c1", 105, 106, 3)]]
    run(tmp_path, REGION, five, ug_plus=ug, name="c", extra=["--nanopare_compat"])
    with pytest.raises(SystemExit):
        run(tmp_path, REGION, five, ug_plus=ug, name="d")
    assert "uG exceeds 5P" in capsys.readouterr().err


def test_chromosome_missing_from_lengths_errors(tmp_path: Path, capsys) -> None:
    with pytest.raises(SystemExit):
        run(tmp_path, [("c9", 0, 10, "g.1", 2, "+")], [[]])
    assert "not in --lengths" in capsys.readouterr().err


def test_replicate_lists_must_match(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        cap_class.main(
            [
                "--ends_plus",
                "a",
                "b",
                "--ends_minus",
                "a",
                "--ug_plus",
                "a",
                "--ug_minus",
                "a",
                "--features",
                "f",
                "--lengths",
                "l",
                "--prefix",
                "p",
            ]
        )


def test_literal_chromosome_names(tmp_path: Path) -> None:
    out = run(
        tmp_path,
        [("Su(var)3-9", 0, 10, "g.1", 2, "+")],
        [[("Su(var)3-9", 2, 3, 4)]],
        lengths=(("Su(var)3-9", 50),),
    )
    assert out["noncapped.bed"][0][0] == "Su(var)3-9"


def test_rule11_order_independent(tmp_path: Path) -> None:
    feats = [
        ("c1", 100, 120, "g.1", 2, "+"),
        ("c1", 300, 330, "g.2", 2, "-"),
        ("c1", 500, 520, "g.3", 2, "+"),
    ]
    reps = [[("c1", 105, 106, 2), ("c1", 510, 512, 1)], [("c1", 106, 107, 3)]]
    mins = [[("c1", 305, 306, 2)], [("c1", 320, 321, 1)]]
    ugp = [[("c1", 105, 106, 1)], []]
    a = run(tmp_path, feats, reps, minus=mins, ug_plus=ugp, name="a")
    rng = random.Random(3)
    feats2 = feats[:]
    rng.shuffle(feats2)
    reps2 = [list(reversed(r)) for r in reversed(reps)]
    b = run(
        tmp_path, feats2, reps2, minus=list(reversed(mins)), ug_plus=list(reversed(ugp)), name="b"
    )
    for k in ("plus.bedgraph", "minus.bedgraph", "peaks.bed", "capped.bed", "noncapped.bed"):
        if k == "peaks.bed":  # names follow the input regions
            assert sorted(a[k]) == sorted(b[k])
        else:
            assert a[k] == b[k], k
