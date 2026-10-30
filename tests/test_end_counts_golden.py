"""Golden test for specs/end-counts.md: `--nanopare_compat` against the original's outputs.

Fixture: tests/fixtures/end_counts_golden/ (README there). Every file must match byte for byte.
"""

from pathlib import Path

import pytest

from degrseq import end_counts

FIXTURE = Path(__file__).parent / "fixtures" / "end_counts_golden"

# golden file -> our output suffix
FILES = {
    "5P": {
        "flower_1.5P.plus.bedgraph": "ends.plus.bedgraph",
        "flower_1.5P.minus.bedgraph": "ends.minus.bedgraph",
        "flower_1.uG.plus.bedgraph": "cap.plus.bedgraph",
        "flower_1.uG.minus.bedgraph": "cap.minus.bedgraph",
    },
    "BODY": {
        "flower_1.BODY.plus.bedgraph": "ends.plus.bedgraph",
        "flower_1.BODY.minus.bedgraph": "ends.minus.bedgraph",
        "flower_1.BODY.cov.bedgraph": "cov.bedgraph",
    },
}


@pytest.mark.parametrize("chunk", [end_counts.COVERAGE_CHUNK, 7])
@pytest.mark.parametrize("library", ["5P", "BODY"])
def test_compat_matches_original(
    tmp_path: Path, library: str, chunk: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A small chunk checks that chunked coverage adds keep the original's summation order.
    monkeypatch.setattr(end_counts, "COVERAGE_CHUNK", chunk)
    prefix = tmp_path / library
    end_counts.main(
        [
            "--bam", str(FIXTURE / "input" / f"flower_1.{library}.bam"),
            "--fasta", str(FIXTURE / "input" / "genome.fasta"),
            "--library", library,
            "--prefix", str(prefix),
            "--nanopare_compat",
            "--write_coverage",
        ]
    )  # fmt: skip
    for golden, ours in FILES[library].items():
        expected = (FIXTURE / "golden" / golden).read_bytes()
        assert Path(f"{prefix}.{ours}").read_bytes() == expected, golden
