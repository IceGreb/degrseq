"""Deterministic SAM subsample: all NH > 1 reads, plus unique reads with name hash < fraction."""

import hashlib
import sys
from itertools import groupby


def keep(name: str, records: list[str], fraction: float) -> bool:
    if any("\tNH:i:" in r and int(r.split("\tNH:i:")[1].split("\t")[0]) > 1 for r in records):
        return True
    h = int.from_bytes(hashlib.md5(name.encode()).digest()[:8], "big")
    return h / 2**64 < fraction


def main() -> None:
    path, fraction = sys.argv[1], float(sys.argv[2])
    out = sys.stdout
    with open(path) as fh:
        lines = iter(fh)
        body = []
        for line in lines:
            if line.startswith("@"):
                out.write(line)
            else:
                body.append(line)
                break
        body_iter = (line for chunk in (body, lines) for line in chunk)
        for name, grp in groupby(body_iter, key=lambda r: r.split("\t", 1)[0]):
            records = list(grp)
            if keep(name, records, fraction):
                out.writelines(records)


if __name__ == "__main__":
    main()
