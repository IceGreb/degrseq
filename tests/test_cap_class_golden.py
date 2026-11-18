"""specs/cap-class.md Tests: compat vs golden EndClass outputs; default vs the quirk table."""

import gzip
import json
import shutil
from collections import defaultdict
from pathlib import Path

from degrseq import cap_class

FIXTURES = Path(__file__).parent / "fixtures"
FIVE = FIXTURES / "scale_factor_test_set" / "input"
GOLDEN = FIXTURES / "cap_class_golden"
LENGTHS = FIXTURES / "kde_peaks_golden" / "length.table"
SAMPLES = ("flower_1", "flower_2", "flower_3")


def unzip(src: Path, dst: Path) -> Path:
    with gzip.open(src, "rb") as a, open(dst, "wb") as b:
        shutil.copyfileobj(a, b)
    return dst


def golden(name: str) -> str:
    with gzip.open(GOLDEN / "golden" / f"flower.{name}.gz", "rt") as fh:
        return fh.read()


def per_nt(path: Path) -> dict[tuple[str, int], float]:
    d: dict[tuple[str, int], float] = defaultdict(float)
    for line in path.read_text().splitlines():
        c, a, b, v = line.split("\t")[:4]
        for i in range(int(a), int(b)):
            d[(c, i)] += float(v)
    return d


def default_track(ug: Path, five: Path, strand: str, out: Path) -> Path:
    """Original uG track moved to the 5' end and capped at 5P per nt (one count per read)."""
    shift = 1 if strand == "plus" else -1
    f = per_nt(five)
    rows = []
    for (c, i), v in sorted(per_nt(ug).items()):
        j = i + shift
        w = min(v, f.get((c, j), 0.0))
        if j >= 0 and w:
            rows.append(f"{c}\t{j}\t{j + 1}\t{w}\n")
    out.write_text("".join(rows))
    return out


def run(tmp_path: Path, compat: bool) -> Path:
    args = ["--features", str(GOLDEN / "input" / "flower.rep.bed"), "--lengths", str(LENGTHS)]
    for strand in ("plus", "minus"):
        fives = [
            unzip(FIVE / f"{s}.5P.{strand}_mask.bedgraph.gz", tmp_path / f"{s}.5P.{strand}")
            for s in SAMPLES
        ]
        ugs = [
            unzip(
                GOLDEN / "input" / f"{s}.uG.{strand}_mask.bedgraph.gz",
                tmp_path / f"{s}.uG.{strand}",
            )
            for s in SAMPLES
        ]
        if not compat:
            ugs = [
                default_track(u, f, strand, Path(f"{u}.default"))
                for u, f in zip(ugs, fives, strict=True)
            ]
        args += [f"--ends_{strand}", *map(str, fives), f"--ug_{strand}", *map(str, ugs)]
    prefix = tmp_path / ("compat" if compat else "default")
    cap_class.main(args + ["--prefix", str(prefix)] + (["--nanopare_compat"] if compat else []))
    return prefix


def columns(text: str, cols: tuple[int, ...]) -> list[tuple[str, ...]]:
    return sorted(tuple(r.split("\t")[i] for i in cols) for r in text.splitlines())


def test_compat_matches_golden(tmp_path: Path) -> None:
    prefix = run(tmp_path, compat=True)
    for name in ("plus.bedgraph", "minus.bedgraph", "uG.plus.bedgraph", "uG.minus.bedgraph"):
        assert Path(f"{prefix}.{name}").read_text() == golden(name), name
    assert columns(Path(f"{prefix}.peaks.bed").read_text(), (0, 1, 2, 3, 4, 5)) == columns(
        golden("peaks.bed"), (0, 1, 2, 3, 4, 5)
    )
    for label, n in (("capped", 323), ("noncapped", 279)):
        ours = Path(f"{prefix}.{label}.bed").read_text()
        assert columns(ours, (0, 1, 2, 5, 9)) == columns(golden(f"{label}.bed"), (0, 1, 2, 5, 9))
        assert len(ours.splitlines()) == n


def test_default_vs_golden(tmp_path: Path) -> None:
    prefix = run(tmp_path, compat=False)
    # G2: unrounded sums move no trimmed span or peak
    assert columns(Path(f"{prefix}.peaks.bed").read_text(), (0, 1, 2, 4, 5)) == columns(
        golden("peaks.bed"), (0, 1, 2, 4, 5)
    )
    gold = {}
    for label in ("capped", "noncapped"):
        for r in golden(f"{label}.bed").splitlines():
            f = r.split("\t")
            gold[(f[0], int(f[1]), int(f[2]) - 1, f[5])] = (
                label,
                int(f[4]) - 1,
            )  # quirk 15: end, peak + 1
    ours = {}
    for label in ("capped", "noncapped"):
        for r in Path(f"{prefix}.{label}.bed").read_text().splitlines():
            f = r.split("\t")
            ours[(f[0], int(f[1]), int(f[2]), f[5])] = (label, int(f[4]))
    assert ours.keys() == gold.keys()  # same features; ends 1 nt shorter than golden
    assert all(ours[k][1] == gold[k][1] for k in ours)  # peak offsets: golden - 1
    gained = sorted(k for k in ours if ours[k][0] == "capped" and gold[k][0] == "noncapped")
    lost = sorted(k for k in ours if ours[k][0] == "noncapped" and gold[k][0] == "capped")
    assert len(gained) == 7 and {k[3] for k in gained} == {
        "+"
    }  # + strand cap G now inside the window
    assert lost == [
        ("Ath_chr2", 225715, 225750, "-"),  # golden fraction relied on several G's per read
        (
            "Ath_chr2",
            489670,
            489672,
            "+",
        ),  # its only G belongs to a read starting outside the feature
        ("Ath_chr2", 722622, 722661, "+"),  # several G's per read
    ]
    stats = json.loads(Path(f"{prefix}.cap-class.stats.json").read_text())
    assert stats["capped"] == {"+": 157, "-": 170} and stats["mode"] == "default"
