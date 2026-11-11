"""Unit tests for specs/scale-factor.md. Each test names the rule or edge case."""

import json
import random
from pathlib import Path

import pytest

from degrseq import scale_factor

# One + gene, region chr1:0-100 (len 100)
REGIONS = "chr1\t0\t100\tgA\t0\t+\n"


def run(
    tmp_path: Path,
    regions: str = REGIONS,
    ends_plus: str = "chr1\t5\t6\t10\n",
    ends_minus: str = "",
    body_plus: str = "",
    body_minus: str = "chr1\t50\t51\t20\n",
    extra: list[str] | None = None,
) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    files = {
        "regions": regions,
        "ends_plus": ends_plus,
        "ends_minus": ends_minus,
        "body_plus": body_plus,
        "body_minus": body_minus,
    }
    args = []
    for name, text in files.items():
        path = tmp_path / f"{name}.txt"
        path.write_text(text)
        args += [f"--{name}", str(path)]
    prefix = tmp_path / "out"
    scale_factor.main(args + ["--prefix", str(prefix), "--write_genes"] + (extra or []))
    return {
        "factor": Path(f"{prefix}.scale_factor.txt").read_text(),
        "stats": json.loads(Path(f"{prefix}.scale-factor.stats.json").read_text()),
        "genes": Path(f"{prefix}.scale-factor.genes.tsv").read_text(),
    }


def genes(result: dict) -> dict[str, tuple[int, float, float]]:
    rows = [line.split("\t") for line in result["genes"].splitlines()[1:]]
    return {r[0]: (int(r[1]), float(r[2]), float(r[3])) for r in rows}


# Rules 6-7: E=10, B=20, L̄=100 → 2·200/100 · 20/10 = 8
def test_formula(tmp_path: Path) -> None:
    result = run(tmp_path)
    assert result["factor"] == "8.0000\n"
    s = result["stats"]
    assert (s["end_total"], s["body_total"], s["mean_gene_length"]) == (10, 20, 100)
    assert s["raw_factor"] == 8 and not s["capped"]


# Rule 1: reverse (default) pairs + genes with body_minus; forward with body_plus
def test_body_strand(tmp_path: Path) -> None:
    kwargs = {"body_plus": "chr1\t50\t51\t40\n", "body_minus": "chr1\t50\t51\t20\n"}
    assert run(tmp_path / "r", **kwargs)["stats"]["body_total"] == 20
    fwd = run(tmp_path / "f", **kwargs, extra=["--body_strand", "forward"])
    assert fwd["stats"]["body_total"] == 40


def test_minus_gene_strand_matching(tmp_path: Path) -> None:
    result = run(
        tmp_path,
        regions="chr1\t0\t100\tgM\t0\t-\n",
        ends_plus="chr1\t5\t6\t99\n",
        ends_minus="chr1\t5\t6\t10\n",
        body_plus="chr1\t50\t51\t20\n",
        body_minus="chr1\t50\t51\t77\n",
    )
    assert (result["stats"]["end_total"], result["stats"]["body_total"]) == (10, 20)


# Rules 2-4: a position shared by 2 genes counts once in E/B, split in half per gene
def test_shared_position_two_genes(tmp_path: Path) -> None:
    regions = REGIONS + "chr1\t40\t200\tgB\t0\t+\n"
    result = run(tmp_path, regions=regions)
    assert result["stats"]["body_total"] == 20
    g = genes(result)
    assert g["gA"][2] == 10 and g["gB"][2] == 10


def test_shared_position_three_genes(tmp_path: Path) -> None:
    regions = REGIONS + "chr1\t40\t200\tgB\t0\t+\nchr1\t50\t60\tgC\t0\t+\n"
    result = run(tmp_path, regions=regions, body_minus="chr1\t50\t51\t30\n")
    assert result["stats"]["body_total"] == 30
    assert {k: v[2] for k, v in genes(result).items()} == {"gA": 10, "gB": 10, "gC": 10}


# Rule 2: opposite-strand gene does not own the position
def test_ownership_is_stranded(tmp_path: Path) -> None:
    regions = REGIONS + "chr1\t0\t100\tgM\t0\t-\n"
    g = genes(run(tmp_path, regions=regions))
    assert g["gA"][2] == 20 and g["gM"][2] == 0


# Rule 5: length = positions in the region; several rows per gene
def test_gene_length_multiple_rows(tmp_path: Path) -> None:
    regions = "chr1\t0\t100\tgA\t0\t+\nchr1\t300\t350\tgA\t0\t+\n"
    result = run(tmp_path, regions=regions)
    assert genes(result)["gA"][0] == 150
    assert result["stats"]["mean_gene_length"] == 150


# Rule 6: L̄ weighted by BODY counts, genes with b = 0 excluded
def test_mean_length_weighting(tmp_path: Path) -> None:
    regions = REGIONS + "chr1\t1000\t1400\tgB\t0\t+\nchr1\t5000\t6000\tgZ\t0\t+\n"
    body = "chr1\t50\t51\t20\nchr1\t1100\t1101\t60\n"
    result = run(tmp_path, regions=regions, body_minus=body)
    # 80 / (20/100 + 60/400)
    assert result["stats"]["mean_gene_length"] == pytest.approx(80 / 0.35)
    assert result["stats"]["genes_with_body"] == 2
    assert result["stats"]["genes_in_regions"] == 3


# Rule 3: values outside regions ignored, reported in stats
def test_outside_values(tmp_path: Path) -> None:
    result = run(
        tmp_path,
        ends_plus="chr1\t5\t6\t10\nchr1\t500\t501\t7\nchrX\t1\t2\t3\n",
        body_minus="chr1\t50\t51\t20\nchr1\t98\t102\t1\n",
    )
    s = result["stats"]
    assert s["end_total"] == 10 and s["body_total"] == 22
    assert s["values_outside_regions"] == {"end": 10, "body": 2}
    assert s["chromosomes_without_regions"] == ["chrX"]


# Inputs: a row spanning n nt counts at each nt; merged == per-nt
def test_merged_rows_equal_per_nt(tmp_path: Path) -> None:
    merged = run(tmp_path / "m", body_minus="chr1\t50\t55\t0.1\n")
    per_nt = run(
        tmp_path / "p", body_minus="".join(f"chr1\t{i}\t{i + 1}\t0.1\n" for i in range(50, 55))
    )
    assert merged["stats"] == per_nt["stats"]
    assert merged["genes"] == per_nt["genes"]


# Rule 8: cap and rounding
def test_cap(tmp_path: Path) -> None:
    result = run(tmp_path, extra=["--max_scale_factor", "5"])
    assert result["factor"] == "5.0000\n"
    assert result["stats"]["capped"] and result["stats"]["raw_factor"] == 8


def test_four_decimals(tmp_path: Path) -> None:
    # raw = 4 · 20/30 = 2.6667
    assert run(tmp_path, ends_plus="chr1\t5\t6\t30\n")["factor"] == "2.6667\n"


# Rule 9: row order does not change any output byte
def test_row_order_invariant(tmp_path: Path) -> None:
    rng = random.Random(1)
    regions = [f"chr1\t{i * 50}\t{i * 50 + 80}\tg{i}\t0\t+\n" for i in range(20)]
    body = [f"chr1\t{i}\t{i + 1}\t{rng.random():.2f}\n" for i in range(0, 1000, 7)]
    ends = [f"chr1\t{i}\t{i + 1}\t{rng.random():.2f}\n" for i in range(0, 1000, 3)]
    a = run(tmp_path / "a", "".join(regions), "".join(ends), "", "", "".join(body))
    for lines in (regions, body, ends):
        rng.shuffle(lines)
    b = run(tmp_path / "b", "".join(regions), "".join(ends), "", "", "".join(body))
    assert a == b


def test_genes_sorted_bytewise(tmp_path: Path) -> None:
    regions = "chr1\t0\t100\tgb\t0\t+\nchr1\t200\t300\tgA\t0\t+\nchr1\t400\t500\tga\t0\t+\n"
    rows = run(tmp_path, regions=regions)["genes"].splitlines()
    assert [r.split("\t")[0] for r in rows[1:]] == ["gA", "ga", "gb"]


# Edge cases
@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"ends_plus": ""}, "E = 0"),
        ({"body_minus": ""}, "B = 0"),
        ({"body_minus": "chr1\t50\t51\t-1\n"}, "finite"),
        ({"body_minus": "chr1\t50\t51\tnan\n"}, "finite"),
        ({"regions": REGIONS + "chr1\t90\t120\tgA\t0\t+\n"}, "overlap"),
        ({"regions": "chr1\t0\t100\tgA\t0\t.\n"}, "strand"),
        ({"regions": "chr1\t0\t100\tgA\n"}, "BED6"),
        ({"ends_plus": "chr1\t6\t5\t1\n"}, "start"),
    ],
)
def test_errors(tmp_path: Path, capsys: pytest.CaptureFixture, kwargs: dict, message: str) -> None:
    with pytest.raises(SystemExit) as exit_info:
        run(tmp_path, **kwargs)
    assert message in str(exit_info.value)


def test_overlapping_rows_of_different_genes_allowed(tmp_path: Path) -> None:
    run(tmp_path, regions=REGIONS + "chr1\t90\t120\tgB\t0\t+\n")


def test_abutting_rows_of_one_gene(tmp_path: Path) -> None:
    regions = "chr1\t0\t50\tgA\t0\t+\nchr1\t50\t100\tgA\t0\t+\n"
    result = run(tmp_path, regions=regions)
    assert result["factor"] == "8.0000\n"
    assert genes(result)["gA"][0] == 100
