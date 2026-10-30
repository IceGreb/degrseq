"""Golden test: `cleavage-stats --nanopare_compat` vs the original's step-2 output (spec Tests).

adj.p.val may differ from R in the 15th digit (scipy's chi-square tail vs R's pchisq, DL-052);
every other field must match as text.
"""

from pathlib import Path

import pytest

from degrseq import cleavage_stats

FIXTURE = Path(__file__).parent / "fixtures" / "cleavage_stats_golden"
INPUT = FIXTURE / "input"


def run(tmp_path: Path, sample: str, extra: list[str]) -> Path:
    prefix = tmp_path / sample
    cleavage_stats.main(
        [
            "--real_sites_w20", str(INPUT / f"{sample}.real.w20.tsv"),
            "--real_sites_w50", str(INPUT / f"{sample}.real.w50.tsv"),
            "--shuffled_sites_w20", str(INPUT / f"{sample}.shuffled.w20.tsv.gz"),
            "--shuffled_sites_w50", str(INPUT / f"{sample}.shuffled.w50.tsv.gz"),
            "--srna_table", str(INPUT / "srna_table.tsv"),
            "--precursors", str(INPUT / "precursors.tsv"),
            "--validated_targets", str(INPUT / "validated_targets.tsv"),
            "--gene_aliases", str(INPUT / "gene_aliases.tsv"),
            "--prefix", str(prefix),
            *extra,
        ]
    )  # fmt: skip
    return prefix


@pytest.mark.parametrize("sample", ["flower_1", "flower_2", "flower_3"])
@pytest.mark.parametrize("cls", ["miRNA", "tasiRNA"])
def test_compat_matches_original(tmp_path: Path, sample: str, cls: str) -> None:
    prefix = run(tmp_path, sample, ["--nanopare_compat"])
    ours = rows(Path(f"{prefix}.{cls}.detected_sites.tsv"))
    golden = rows(FIXTURE / "golden" / f"{sample}.{cls}.detected_sites.tsv")
    assert [r[:-1] for r in ours] == [r[:-1] for r in golden]
    for mine, theirs in zip(ours[1:], golden[1:], strict=True):
        assert float(mine[-1]) == pytest.approx(float(theirs[-1]), rel=1e-12, abs=0)


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()]
