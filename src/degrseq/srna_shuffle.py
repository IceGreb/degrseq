"""Shuffled small-RNA sets: the null model for cleavage calling.

Implements specs/srna-shuffle.md (approved 30/10/2026).
"""

import argparse
import json
import random
import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TextIO

SCHEMA_VERSION = 1
EXPECTED_CHARS = set("ACGTUNacgtun")


class SrnaShuffleError(Exception):
    pass


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="degrseq srna-shuffle",
        description="Shuffle each sRNA N times, keeping length and composition.",
    )
    parser.add_argument("--fasta", required=True, help="sRNA FASTA")
    parser.add_argument("--n_shuffles", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--layout", choices=("combined", "per_set"), default="combined")
    parser.add_argument("--prefix", help="output file prefix (default: input file stem)")
    args = parser.parse_args(argv)
    if args.n_shuffles < 1:
        parser.error("--n_shuffles must be >= 1")
    if args.prefix is None:
        args.prefix = Path(args.fasta).stem
    return args


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Rule 1 and edge cases: (ID, sequence) in input order; multi-line records joined."""
    records: list[tuple[str, list[str]]] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise SrnaShuffleError("header without an ID")
                records.append((fields[0], []))
            elif not records:
                raise SrnaShuffleError("sequence before the first header")
            else:
                records[-1][1].append(line)
    seen: set[str] = set()
    result = []
    for name, parts in records:
        if name in seen:
            raise SrnaShuffleError(f"duplicate ID {name}")
        seen.add(name)
        seq = "".join(parts)
        if not seq:
            raise SrnaShuffleError(f"record {name} has an empty sequence")
        result.append((name, seq))
    return result


def shuffled_sets(
    records: list[tuple[str, str]], n_shuffles: int, seed: int
) -> Iterator[tuple[int, list[tuple[str, str, str]]]]:
    """Rules 2, 5: one random stream, set-major; yields (j, [(ID, input, shuffle)])."""
    rng = random.Random(seed)
    for j in range(n_shuffles):
        batch = []
        for name, seq in records:
            chars = list(seq)
            rng.shuffle(chars)
            batch.append((name, seq, "".join(chars)))
        yield j, batch


def write_set(out: TextIO, j: int, batch: list[tuple[str, str, str]]) -> int:
    """Rule 4; returns the number of shuffles identical to their input (rule 6)."""
    for name, _, shuffled in batch:
        out.write(f">{name}_shuffled_{j}\n{shuffled}\n")
    return sum(shuffled == seq for _, seq, shuffled in batch)


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    try:
        records = read_fasta(args.fasta)
    except (OSError, SrnaShuffleError) as error:
        sys.exit(f"degrseq srna-shuffle: error: {error}")
    for name, seq in records:
        if set(seq) - EXPECTED_CHARS:
            print(
                f"degrseq srna-shuffle: warning: {name} has characters outside ACGTUN; "
                "shuffled as-is",
                file=sys.stderr,
            )

    n_identical = 0
    sets = shuffled_sets(records, args.n_shuffles, args.seed)
    if args.layout == "combined":
        with open(f"{args.prefix}.shuffled.fa", "w") as out:
            for j, batch in sets:
                n_identical += write_set(out, j, batch)
    else:
        for j, batch in sets:
            with open(f"{args.prefix}.{j}.shuffled.fa", "w") as out:
                n_identical += write_set(out, j, batch)

    stats = {
        "schema_version": SCHEMA_VERSION,
        "n_srnas": len(records),
        "n_shuffles": args.n_shuffles,
        "seed": args.seed,
        "n_identical_to_input": n_identical,
    }
    stats_path = Path(f"{args.prefix}.srna-shuffle.stats.json")
    stats_path.write_text(json.dumps(stats, indent=2) + "\n")
