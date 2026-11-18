"""specs/kde-peaks.md Tests: compat vs golden EndGraph outputs; default vs the quirk table."""

import gzip
import json
import shutil
from pathlib import Path

import pytest

from degrseq import kde_peaks

FIXTURES = Path(__file__).parent / "fixtures"
INPUT = FIXTURES / "scale_factor_test_set" / "input"
GOLDEN = FIXTURES / "kde_peaks_golden"
SAMPLES = ("flower_1", "flower_2", "flower_3")
# features called / kept (plus + minus), default mode; quirk table + G3 (spec)
DEFAULT_COUNTS = {"flower_1": (1090, 982), "flower_2": (1283, 1126), "flower_3": (1219, 1054)}


def unzip(src: Path, dst: Path) -> str:
    with gzip.open(src, "rb") as a, open(dst, "wb") as b:
        shutil.copyfileobj(a, b)
    return str(dst)


def run(sample: str, tmp_path: Path, compat: bool) -> Path:
    factors = dict(
        line.split("\t") for line in (GOLDEN / "scale_factors.tsv").read_text().splitlines()
    )
    sf = tmp_path / "sf.txt"
    sf.write_text(factors[sample] + "\n")
    args = ["--scale_factor", str(sf), "--lengths", str(GOLDEN / "length.table"), "--write_density"]
    for flag, name in [
        ("--ends_plus", "5P.plus_mask"),
        ("--ends_minus", "5P.minus_mask"),
        ("--body_plus", "BODY.plus"),
        ("--body_minus", "BODY.minus"),
    ]:
        args += [flag, unzip(INPUT / f"{sample}.{name}.bedgraph.gz", tmp_path / f"{name}.bedgraph")]
    prefix = tmp_path / ("compat" if compat else "default")
    kde_peaks.main(args + ["--prefix", str(prefix)] + (["--nanopare_compat"] if compat else []))
    return prefix


def sorted_rows(text: str) -> list[list[str]]:
    rows = [line.split("\t") for line in text.splitlines()]
    return sorted(rows, key=lambda r: (r[0], int(r[1]), r[5]))


@pytest.mark.parametrize("sample", SAMPLES)
def test_compat_matches_golden(sample: str, tmp_path: Path) -> None:
    prefix = run(sample, tmp_path, compat=True)
    for label in ("plus", "minus"):
        with gzip.open(GOLDEN / "golden" / f"{sample}.5P.{label}_smooth.bedgraph.gz", "rt") as fh:
            golden = fh.read()
        assert Path(f"{prefix}.density_{label}.bedgraph").read_text() == golden
    with gzip.open(GOLDEN / "golden" / f"{sample}.end_features.bed.gz", "rt") as fh:
        golden_bed = fh.read()
    assert sorted_rows(Path(f"{prefix}.end_features.bed").read_text()) == sorted_rows(golden_bed)


@pytest.mark.parametrize("sample", SAMPLES)
def test_default_counts(sample: str, tmp_path: Path) -> None:
    prefix = run(sample, tmp_path, compat=False)
    stats = json.loads(Path(f"{prefix}.kde-peaks.stats.json").read_text())
    called = stats["plus"]["features_called"] + stats["minus"]["features_called"]
    kept = stats["plus"]["features_kept"] + stats["minus"]["features_kept"]
    assert (called, kept) == DEFAULT_COUNTS[sample]
    rows = Path(f"{prefix}.end_features.bed").read_text().splitlines()
    assert len(rows) == kept
    # quirk 12 fixed: BED end is exclusive, so every feature spans >= min_length nt
    assert all(int(r.split("\t")[2]) - int(r.split("\t")[1]) >= 10 for r in rows)
