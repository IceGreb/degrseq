"""Toy genome, SAM builder and runner for the end-counts unit tests."""

import json
import random
from pathlib import Path

import pysam

from degrseq import end_counts

CHROM_LENGTH = 300


def _toy_genome() -> dict[str, str]:
    # No G in the genome, so an untemplated G is never reclaimed.
    rng = random.Random(1)
    chr1 = "".join(rng.choice("ACT") for _ in range(CHROM_LENGTH))
    chr2 = list("".join(rng.choice("ACT") for _ in range(CHROM_LENGTH)))
    chr2[100] = "N"
    return {"chr1": chr1, "chr2": "".join(chr2)}


GENOME = _toy_genome()


def sam(
    name: str,
    flag: int,
    chrom: str,
    pos: int,
    cigar: str,
    seq: str | None = None,
    tags: str = "",
) -> str:
    """One SAM line; `pos` is 0-based. `seq` defaults to A's of the query length."""
    if seq is None:
        n, length = "", 0
        for char in cigar:
            if char.isdigit():
                n += char
            else:
                if char in "MIS=X":
                    length += int(n)
                n = ""
        seq = "A" * length
    fields = [name, str(flag), chrom, str(pos + 1), "255", cigar, "*", "0", "0", seq, "*"]
    line = "\t".join(fields)
    return line + ("\t" + tags if tags else "")


def build_inputs(tmp_path: Path, records: list[str], sort_order: str = "unsorted") -> tuple:
    tmp_path.mkdir(parents=True, exist_ok=True)
    fasta = tmp_path / "genome.fa"
    fasta.write_text("".join(f">{c}\n{s}\n" for c, s in GENOME.items()))
    pysam.faidx(str(fasta))
    header = f"@HD\tVN:1.6\tSO:{sort_order}\n" + "".join(
        f"@SQ\tSN:{c}\tLN:{len(s)}\n" for c, s in GENOME.items()
    )
    sam_path = tmp_path / "reads.sam"
    sam_path.write_text(header + "".join(r + "\n" for r in records))
    bam_path = tmp_path / "reads.bam"
    with (
        pysam.AlignmentFile(str(sam_path), "r") as sam_in,
        pysam.AlignmentFile(str(bam_path), "wb", template=sam_in) as bam_out,
    ):
        for record in sam_in:
            bam_out.write(record)
    return bam_path, fasta


def read_bedgraph(path: Path) -> dict[tuple[str, int], str]:
    rows = {}
    for line in path.read_text().splitlines():
        chrom, start, end, value = line.split("\t")
        rows[(chrom, int(start))] = value
        assert int(end) == int(start) + 1
    return rows


def run(
    tmp_path: Path, records: list[str], library: str = "5P", extra: list[str] | None = None
) -> dict:
    bam, fasta = build_inputs(tmp_path, records)
    prefix = tmp_path / "out"
    end_counts.main(
        ["--bam", str(bam), "--fasta", str(fasta), "--library", library, "--prefix", str(prefix)]
        + (extra or [])
    )
    result = {"stats": json.loads(Path(f"{prefix}.end-counts.stats.json").read_text())}
    tracks = ["ends", "cap"] if library == "5P" else ["ends"]
    for track in tracks:
        for label, strand in (("plus", "+"), ("minus", "-")):
            result[(track, strand)] = read_bedgraph(Path(f"{prefix}.{track}.{label}.bedgraph"))
    result["prefix"] = prefix
    return result


def fates(result: dict) -> dict[str, int]:
    return {k: v for k, v in result["stats"]["read_fates"].items() if v}
