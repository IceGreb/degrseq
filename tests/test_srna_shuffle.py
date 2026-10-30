"""Unit tests for specs/srna-shuffle.md. Each test names the rule or edge case."""

import json
import os
import subprocess
import sys
from collections import Counter
from itertools import permutations
from pathlib import Path

import pytest
from scipy.stats import chisquare

from degrseq import srna_shuffle

FASTA = ">miR1\nUCGGACCAGGCUUCAUUCCCC\n>miR2\nTGAGCCGAACCAATATCACTC\n>tas3\nttcttgaccttgtaagacccc\n"


def run(tmp_path: Path, fasta: str = FASTA, extra: list[str] | None = None) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "srnas.fa"
    path.write_text(fasta)
    prefix = tmp_path / "out"
    args = ["--fasta", str(path), "--prefix", str(prefix), "--n_shuffles", "5"]
    srna_shuffle.main(args + (extra or []))
    stats = json.loads(Path(f"{prefix}.srna-shuffle.stats.json").read_text())
    return {"prefix": prefix, "stats": stats}


def records(text: str) -> list[tuple[str, str]]:
    lines = text.splitlines()
    return list(zip(lines[0::2], lines[1::2], strict=True))


def combined(result: dict) -> str:
    return Path(f"{result['prefix']}.shuffled.fa").read_text()


# Rules 2, 4
def test_headers_set_major_one_line(tmp_path: Path) -> None:
    out = records(combined(run(tmp_path)))
    expected = [f">{n}_shuffled_{j}" for j in range(5) for n in ("miR1", "miR2", "tas3")]
    assert [h for h, _ in out] == expected


# Rule 3
def test_length_and_composition_and_case_kept(tmp_path: Path) -> None:
    inputs = dict(records(FASTA))
    for header, seq in records(combined(run(tmp_path))):
        original = inputs[">" + header[1:].rsplit("_shuffled_", 1)[0]]
        assert Counter(seq) == Counter(original)


# Rule 5
def test_all_orders_equally_likely() -> None:
    orders = Counter(
        batch[0][2] for _, batch in srna_shuffle.shuffled_sets([("s", "ACGT")], 24_000, seed=0)
    )
    counts = [orders["".join(p)] for p in permutations("ACGT")]
    assert sum(counts) == 24_000
    assert chisquare(counts).pvalue > 0.001


def test_same_seed_same_bytes(tmp_path: Path) -> None:
    first = combined(run(tmp_path / "a", extra=["--seed", "7"]))
    second = combined(run(tmp_path / "b", extra=["--seed", "7"]))
    assert first == second


def test_different_seeds_differ(tmp_path: Path) -> None:
    assert combined(run(tmp_path / "a")) != combined(run(tmp_path / "b", extra=["--seed", "1"]))


def test_per_set_layout_holds_the_same_sequences(tmp_path: Path) -> None:
    whole = combined(run(tmp_path / "c"))
    result = run(tmp_path / "p", extra=["--layout", "per_set"])
    parts = [Path(f"{result['prefix']}.{j}.shuffled.fa").read_text() for j in range(5)]
    assert "".join(parts) == whole
    assert not Path(f"{result['prefix']}.shuffled.fa").exists()


def test_output_independent_of_hash_seed(tmp_path: Path) -> None:
    path = tmp_path / "srnas.fa"
    path.write_text(FASTA)
    outputs = []
    for hash_seed in ("0", "123"):
        prefix = tmp_path / f"h{hash_seed}"
        env = {**os.environ, "PYTHONHASHSEED": hash_seed}
        subprocess.run(
            [sys.executable, "-m", "degrseq", "srna-shuffle", "--fasta", str(path),
             "--prefix", str(prefix), "--n_shuffles", "3"],
            check=True, env=env,
        )  # fmt: skip
        outputs.append(Path(f"{prefix}.shuffled.fa").read_bytes())
    assert outputs[0] == outputs[1]


# Rule 6
def test_identical_shuffles_kept_and_counted(tmp_path: Path) -> None:
    result = run(tmp_path, ">polyA\nAAAAAAAAAAAAAAAAAAAA\n")
    assert records(combined(result))[0][1] == "A" * 20
    assert result["stats"]["n_identical_to_input"] == 5


# Rule 1, edge cases
def test_multi_line_fasta(tmp_path: Path) -> None:
    out = records(combined(run(tmp_path, ">miR1\nUCGGACCAGG\nCUUCAUUCCCC\n\n")))
    assert len(out[0][1]) == 21


def test_description_dropped_from_header(tmp_path: Path) -> None:
    out = records(combined(run(tmp_path, ">miR1 MIMAT0000001 ath\nUCGGACCAGG\n")))
    assert out[0][0] == ">miR1_shuffled_0"


@pytest.mark.parametrize(
    ("fasta", "message"),
    [
        (">miR1\nACGT\n>miR1\nACGA\n", "duplicate ID miR1"),
        (">miR1\n>miR2\nACGT\n", "miR1 has an empty sequence"),
        ("ACGT\n>miR1\nACGT\n", "before the first header"),
    ],
)
def test_bad_input_exits(tmp_path: Path, fasta: str, message: str) -> None:
    with pytest.raises(SystemExit, match=message):
        run(tmp_path, fasta)


@pytest.mark.parametrize("n", ["0", "-3"])
def test_n_shuffles_must_be_positive(tmp_path: Path, n: str) -> None:
    with pytest.raises(SystemExit):
        srna_shuffle.main(["--fasta", str(tmp_path / "x.fa"), "--n_shuffles", n])


def test_unexpected_characters_warn(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    run(tmp_path, ">odd\nACGTRY\n")
    assert "odd has characters outside ACGTUN" in capsys.readouterr().err


def test_default_prefix_is_input_stem(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "anno.mir.tas.fa"
    path.write_text(FASTA)
    monkeypatch.chdir(tmp_path)
    srna_shuffle.main(["--fasta", str(path), "--n_shuffles", "1"])
    assert (tmp_path / "anno.mir.tas.shuffled.fa").exists()


def test_stats(tmp_path: Path) -> None:
    stats = run(tmp_path, extra=["--seed", "3"])["stats"]
    assert stats == {
        "schema_version": 1,
        "n_srnas": 3,
        "n_shuffles": 5,
        "seed": 3,
        "n_identical_to_input": 0,
    }
