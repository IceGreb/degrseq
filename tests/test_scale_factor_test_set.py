"""specs/scale-factor.md Tests: default mode on nanoPARE's test set vs the audit's values."""

import gzip
import json
import shutil
from pathlib import Path

import pytest

from degrseq import scale_factor

FIXTURE = Path(__file__).parent / "fixtures" / "scale_factor_test_set" / "input"
EXPECTED = {
    "flower_1": ("0.3871", 138729, 160859),
    "flower_2": ("0.1903", 214506, 121465),
    "flower_3": ("0.2103", 219022, 139615),
}


def unzip(name: str, tmp_path: Path) -> str:
    out = tmp_path / name.removesuffix(".gz")
    with gzip.open(FIXTURE / name, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return str(out)


@pytest.mark.parametrize("sample", sorted(EXPECTED))
def test_test_set(sample: str, tmp_path: Path) -> None:
    factor, end_total, body_total = EXPECTED[sample]
    args = ["--regions", unzip("regions.bed.gz", tmp_path), "--prefix", str(tmp_path / sample)]
    for flag, name in [
        ("--ends_plus", "5P.plus_mask"),
        ("--ends_minus", "5P.minus_mask"),
        ("--body_plus", "BODY.plus"),
        ("--body_minus", "BODY.minus"),
    ]:
        args += [flag, unzip(f"{sample}.{name}.bedgraph.gz", tmp_path)]
    scale_factor.main(args)
    assert Path(f"{tmp_path / sample}.scale_factor.txt").read_text() == factor + "\n"
    stats = json.loads(Path(f"{tmp_path / sample}.scale-factor.stats.json").read_text())
    assert round(stats["end_total"]) == end_total
    assert round(stats["body_total"]) == body_total
    assert not stats["capped"]
