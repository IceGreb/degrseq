"""Gene regions for the scale-factor fixture: exons + 100-nt flanks, merged per gene."""

import sys
from collections import defaultdict

BUFFER = 100
EXON_TYPES = ("exon", "pseudogenic_exon")


def main(gff: str) -> None:
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    info: dict[str, tuple[str, str, str]] = {}
    with open(gff, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in EXON_TYPES:
                continue
            for attr in fields[8].split(";"):
                if attr.startswith("Parent="):
                    for tx in attr[7:].split(","):
                        exons[tx].append((int(fields[3]) - 1, int(fields[4])))
                        info[tx] = (fields[0], fields[6], tx.split(".")[0])

    intervals: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    for tx, ex in exons.items():
        ex.sort()
        flanks = [(max(ex[0][0] - BUFFER, 0), ex[0][0]), (ex[-1][1], ex[-1][1] + BUFFER)]
        intervals[info[tx]] += ex + flanks

    for (chrom, strand, gene), ivs in sorted(intervals.items()):
        ivs.sort()
        merged: list[list[int]] = []
        for start, end in ivs:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        for start, end in merged:
            print(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}")


if __name__ == "__main__":
    main(sys.argv[1])
