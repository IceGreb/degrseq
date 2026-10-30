"""Significant sRNA cleavage sites: empirical p-values against shuffled sRNAs.

Implements specs/cleavage-stats.md (approved 30/10/2026): default mode and `--nanopare_compat`.
"""

import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.stats import chi2

SCHEMA_VERSION = 1
WINDOWS = (20, 50)
CLASSES = ("miRNA", "tasiRNA")
NAME_PREFIXES = {"miR": "miRNA", "TAS": "tasiRNA"}
FAMILY_CHARS = {"miRNA": 6, "tasiRNA": 4}
SITE_COLUMNS = ("sRNA", "Transcript", "Slice.site", "AllenScore", "flankMax.rpm", "slice.site.rpm")
HEADER = (
    "target.transcript",
    "sRNA",
    "slice.site",
    "Allen.score",
    "slice.site.rpm",
    "fold.change",
    "verified.",
    "common.name",
    "interaction",
    "adj.p.val",
)
MISSING = ("", "NA")
Null = dict[str, tuple[np.ndarray, np.ndarray]]


class CleavageStatsError(Exception):
    """Input that violates the spec's preconditions."""


@dataclass
class Site:
    srna: str
    transcript: str
    slice_site: int
    allen: float
    rpm: float
    fc: float
    flank: float = 0.0
    n_null: int = 0
    p_fc: float = 1.0
    p_allen: float = 1.0
    p: float = 1.0
    adj: float = 1.0

    @property
    def gene(self) -> str:
        return self.transcript.split(".", 1)[0]


@dataclass
class Annotation:
    classes: dict[str, str]
    families: dict[str, str]
    precursors: dict[str, set[str]]
    validated: list[tuple[str, str]]
    aliases: dict[str, str]
    has_validated: bool
    _validated_keys: dict[str, set[tuple[str, str]]] = field(default_factory=dict)

    def family(self, srna: str, cls: str) -> str:
        """Rule 13."""
        return self.families.get(srna) or srna[: FAMILY_CHARS[cls]]

    def common(self, gene: str) -> str:
        return self.aliases.get(gene, "NA")

    def verified(self, site: Site, cls: str) -> str:
        """Rule 11: a validated entry matches by its name or its family."""
        if cls not in self._validated_keys:
            self._validated_keys[cls] = {
                (gene, key) for s, gene in self.validated for key in (s, self.family(s, cls))
            }
        hit = (site.gene, self.family(site.srna, cls)) in self._validated_keys[cls]
        return "yes" if hit else "no"


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="degrseq cleavage-stats",
        description="Call significant sRNA cleavage sites against shuffled-sRNA nulls.",
    )
    for w in WINDOWS:
        parser.add_argument(f"--real_sites_w{w}", required=True, help=f"real step-1 table, WIN {w}")
    for w in WINDOWS:
        parser.add_argument(
            f"--shuffled_sites_w{w}", required=True, help=f"pooled shuffled step-1 table, WIN {w}"
        )
    parser.add_argument("--srna_table", required=True, help="TSV: sRNA[, class, family]")
    parser.add_argument("--precursors", help="TSV: sRNA, precursor")
    parser.add_argument("--validated_targets", help="TSV: sRNA, targetAGI")
    parser.add_argument("--gene_aliases", help="TSV: locus_name, symbol")
    parser.add_argument("--prefix", required=True, help="output file prefix")
    parser.add_argument("--fc_min", type=float, default=1.0)
    parser.add_argument("--padj_max", type=float, default=0.05)
    parser.add_argument("--rpm_min", type=float, default=0.1)
    parser.add_argument("--collapse_nt", type=int, default=5)
    parser.add_argument("--nanopare_compat", action="store_true")
    args = parser.parse_args(argv)
    if args.collapse_nt < 0:
        parser.error("--collapse_nt must be >= 0")
    return args


def iter_tsv(path: str, required: Sequence[str]) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.endswith(".gz") else open
    try:
        with opener(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            names = [n.strip() for n in reader.fieldnames or []]
            missing = [c for c in required if c not in names]
            if missing:
                raise CleavageStatsError(f"{path}: missing column(s) {', '.join(missing)}")
            reader.fieldnames = names
            for row in reader:
                yield {k: (v or "").strip() for k, v in row.items() if k is not None}
    except (csv.Error, UnicodeDecodeError, EOFError, gzip.BadGzipFile) as error:
        raise CleavageStatsError(f"{path}: {error}") from None


def number(value: str, path: str, column: str, shuffled: bool) -> float | None:
    """Rule 5: missing (incl. NaN) is None in shuffled rows, an error in real rows."""
    try:
        x = float(value) if value not in MISSING else math.nan
    except ValueError:
        raise CleavageStatsError(f"{path}: {column} {value!r} is not a number") from None
    if math.isnan(x):
        if shuffled:
            return None
        raise CleavageStatsError(f"{path}: missing {column} in a real-site row")
    if not 0 <= x < math.inf:
        raise CleavageStatsError(f"{path}: {column} {value!r} must be finite and >= 0")
    return x


def fold_change(rpm: float | None, flank: float | None) -> float:
    """Rule 1; 0 if either value is missing (rule 5)."""
    return 0.0 if rpm is None or flank is None else (rpm + 1) / (flank + 1)


def read_sites(path: str) -> list[Site]:
    """Real sites."""
    sites = []
    for row in iter_tsv(path, SITE_COLUMNS):
        allen, flank, rpm = (number(row[c], path, c, False) for c in SITE_COLUMNS[3:])
        try:
            slice_site = int(row["Slice.site"])
        except ValueError:
            raise CleavageStatsError(
                f"{path}: Slice.site {row['Slice.site']!r} is not an integer"
            ) from None
        sites.append(
            Site(
                row["sRNA"],
                row["Transcript"],
                slice_site,
                allen,
                rpm,
                fold_change(rpm, flank),
                flank,
            )
        )
    return sites


def read_null(path: str) -> Null:
    """Rule 5: pooled shuffled rows as fc and Allen arrays per sRNA, in input order."""
    values: dict[str, tuple[list[float], list[float]]] = defaultdict(lambda: ([], []))
    for row in iter_tsv(path, SITE_COLUMNS):
        allen, flank, rpm = (number(row[c], path, c, True) for c in SITE_COLUMNS[3:])
        fcs, allens = values[row["sRNA"]]
        fcs.append(fold_change(rpm, flank))
        allens.append(allen or 0.0)
    return {k: (np.array(f), np.array(a)) for k, (f, a) in values.items()}


def read_annotation(args: argparse.Namespace) -> tuple[Annotation, list[str]]:
    """Rule 13 and the optional label inputs; returns the unclassified sRNAs too."""
    classes: dict[str, str] = {}
    families: dict[str, str] = {}
    unclassified = []
    for row in iter_tsv(args.srna_table, ("sRNA",)):
        name = row["sRNA"]
        if name in classes or name in unclassified:
            raise CleavageStatsError(f"{args.srna_table}: duplicate sRNA {name}")
        cls = row.get("class", "")
        if cls and cls not in CLASSES:
            raise CleavageStatsError(f"{args.srna_table}: {name} has unknown class {cls!r}")
        cls = cls or NAME_PREFIXES.get(name[:3], "")
        if not cls:
            unclassified.append(name)
            continue
        classes[name] = cls
        if row.get("family"):
            families[name] = row["family"]

    precursors: dict[str, set[str]] = defaultdict(set)
    if args.precursors:
        for r in iter_tsv(args.precursors, ("sRNA", "precursor")):
            precursors[r["sRNA"]].add(r["precursor"])
    validated = []
    if args.validated_targets:
        rows = iter_tsv(args.validated_targets, ("sRNA", "targetAGI"))
        validated = [(r["sRNA"], r["targetAGI"]) for r in rows]
    aliases: dict[str, str] = {}
    if args.gene_aliases:
        for r in iter_tsv(args.gene_aliases, ("locus_name", "symbol")):
            aliases.setdefault(r["locus_name"], r["symbol"])
    annotation = Annotation(
        classes, families, precursors, validated, aliases, bool(args.validated_targets)
    )
    return annotation, unclassified


def collapse(sites: list[Site], nt: int) -> list[Site]:
    """Rule 4."""
    ordered = sorted(sites, key=lambda s: (s.allen, -s.rpm, s.srna))
    kept: dict[str, list[Site]] = {}
    for site in ordered:
        group = kept.setdefault(site.transcript, [])
        if all(abs(site.slice_site - k.slice_site) > nt for k in group):
            group.append(site)
    return [s for group in kept.values() for s in group]


def bh(p: np.ndarray) -> np.ndarray:
    """Rule 8, with R's p.adjust arithmetic."""
    n = p.size
    order = np.argsort(-p, kind="stable")
    ranks = np.arange(n, 0, -1)
    adjusted = np.minimum(1.0, np.minimum.accumulate((n / ranks) * p[order]))
    out = np.empty(n)
    out[order] = adjusted
    return out


def window_pvalues(real: list[Site], null: Null, empty_null: set[str]) -> list[Site]:
    """Rules 6 to 8."""
    sorted_null: Null = {}
    for site in real:
        if site.srna not in null:
            empty_null.add(site.srna)
            continue
        if site.srna not in sorted_null:
            sorted_null[site.srna] = tuple(np.sort(v) for v in null[site.srna])
        fcs, allens = sorted_null[site.srna]
        site.n_null = fcs.size
        site.p_fc = float(1 - np.searchsorted(fcs, site.fc, side="right") / fcs.size)
        site.p_allen = float(np.searchsorted(allens, site.allen, side="right") / allens.size)
        with np.errstate(divide="ignore"):
            statistic = -2 * (np.log(site.p_fc) + np.log(site.p_allen))
        site.p = float(chi2.sf(statistic, 4))
    for site, adj in zip(real, bh(np.array([s.p for s in real])), strict=True):
        site.adj = float(adj)
    return real


def r_number(x: float) -> str:
    """Rule 14: R write.table formatting of a double."""
    if x == 0:
        return "0"
    reference = float(f"{x:.14e}")
    for digits in range(1, 16):
        sci = f"{x:.{digits - 1}e}"
        if float(sci) == reference:
            break
    exponent = int(sci.split("e")[1])
    fixed = f"{x:.{max(0, digits - 1 - exponent)}f}"
    return sci if len(fixed) > len(sci) else fixed


def join_windows(
    by_window: dict[int, list[Site]], cls: str, annotation: Annotation, compat: bool
) -> list[tuple]:
    """Rules 10 to 12."""
    keyed: dict[tuple[str, int], dict[int, Site]] = defaultdict(dict)
    for w, sites in by_window.items():
        for site in sites:
            keyed[(site.transcript, site.slice_site)][w] = site
    rows = []
    for (transcript, slice_site), hits in keyed.items():
        first, second = hits.get(WINDOWS[0]), hits.get(WINDOWS[1])
        chosen = first if first and (second is None or first.adj <= second.adj) else second
        verified = [annotation.verified(s, cls) for s in hits.values()]
        common = annotation.common(chosen.gene)
        if compat:
            verified_out = min(verified)
            common_out = common if len(hits) == 2 else "NA"
            srna = (first or second).srna
        else:
            verified_out = "yes" if "yes" in verified else "no"
            common_out = common
            srna = chosen.srna
        target = common if common != "NA" else transcript
        rows.append(
            (transcript, chosen.srna, slice_site, chosen.allen, chosen.rpm, chosen.fc,
             verified_out, common_out, f"{srna}:{target}", chosen.adj)
        )  # fmt: skip
    if compat:
        rows.sort(key=lambda r: f"{r[0]}_{r[2]}")
        rows.sort(key=lambda r: r[9])
    else:
        rows.sort(key=lambda r: (r[9], -r[5], r[1], r[0], r[2]))
    return rows


def format_row(row: tuple) -> str:
    transcript, srna, slice_site, allen, rpm, fc, verified, common, interaction, adj = row
    fields = [transcript, srna, str(slice_site), r_number(allen), r_number(rpm), r_number(fc),
              verified, common, interaction, r_number(adj)]  # fmt: skip
    return "\t".join(fields) + "\n"


def write_table(path: str, rows: list[tuple]) -> None:
    with open(path, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        out.writelines(format_row(r) for r in rows)


WINDOW_HEADER = ("window", "sRNA", "Transcript", "Slice.site", "AllenScore", "flankMax.rpm",
                 "slice.site.rpm", "fold.change", "n_shuffled", "fc.p.val", "allen.p.val",
                 "fishers.combined.p.val", "fishers.combined.p.val.adj")  # fmt: skip
SHUFFLED_HEADER = ("window", "sRNA", "fold.change", "AllenScore")


def write_window_tests(path: str, by_window: dict[int, list[Site]]) -> None:
    """Rule 15: every tested real site per window, with its p-values."""
    with open(path, "w") as out:
        out.write("\t".join(WINDOW_HEADER) + "\n")
        for w in WINDOWS:
            for s in sorted(by_window[w], key=lambda s: (s.srna, s.transcript, s.slice_site)):
                numbers = (s.allen, s.flank, s.rpm, s.fc)
                tail = (s.p_fc, s.p_allen, s.p, s.adj)
                fields = [str(w), s.srna, s.transcript, str(s.slice_site),
                          *map(r_number, numbers), str(s.n_null), *map(r_number, tail)]  # fmt: skip
                out.write("\t".join(fields) + "\n")


def write_shuffled_values(path: str, rows: list[tuple[int, str, float, float]]) -> None:
    """Rule 15: the shuffled values behind each tested sRNA's p-values."""
    with open(path, "w") as out:
        out.write("\t".join(SHUFFLED_HEADER) + "\n")
        for w, srna, fc, allen in rows:
            out.write(f"{w}\t{srna}\t{r_number(fc)}\t{r_number(allen)}\n")


def significant(row: tuple, args: argparse.Namespace) -> bool:
    """Rule 9."""
    return row[5] > args.fc_min and row[9] < args.padj_max and row[4] > args.rpm_min


def run(args: argparse.Namespace) -> dict:
    annotation, unclassified = read_annotation(args)
    for name in unclassified:
        print(
            f"degrseq cleavage-stats: warning: {name} has no class and no miR/TAS prefix; skipped",
            file=sys.stderr,
        )
    empty_null: set[str] = set()
    not_in_table: set[str] = set()
    by_window: dict[str, dict[int, list[Site]]] = {cls: {} for cls in CLASSES}
    shuffled_rows: dict[str, list[tuple[int, str, float, float]]] = {cls: [] for cls in CLASSES}
    window_stats: dict[str, dict[str, dict]] = {cls: {} for cls in CLASSES}
    for w in WINDOWS:
        real = read_sites(getattr(args, f"real_sites_w{w}"))
        not_in_table |= {s.srna for s in real if s.srna not in annotation.classes}
        null = read_null(getattr(args, f"shuffled_sites_w{w}"))
        for cls in CLASSES:
            in_class = [s for s in real if annotation.classes.get(s.srna) == cls]
            kept = [s for s in in_class if s.gene not in annotation.precursors.get(s.srna, ())]
            collapsed = collapse(kept, args.collapse_nt)
            by_window[cls][w] = window_pvalues(collapsed, null, empty_null)
            for srna in sorted({s.srna for s in collapsed} & null.keys()):
                fcs, allens = null[srna]
                shuffled_rows[cls] += [
                    (w, srna, f, a) for f, a in zip(fcs.tolist(), allens.tolist(), strict=True)
                ]
            window_stats[cls][str(w)] = {
                "n_real_rows": len(in_class),
                "n_after_self_self": len(kept),
                "n_after_collapse": len(collapsed),
                "n_srnas_tested": len({s.srna for s in collapsed}),
                "n_shuffled_rows": sum(
                    v[0].size for k, v in null.items() if annotation.classes.get(k) == cls
                ),
            }
        del null

    class_stats = {}
    for cls in CLASSES:
        rows = join_windows(by_window[cls], cls, annotation, args.nanopare_compat)
        hits = [r for r in rows if significant(r, args)]
        write_table(f"{args.prefix}.{cls}.detected_sites.tsv", rows)
        write_table(f"{args.prefix}.{cls}.significant_sites.tsv", hits)
        write_window_tests(f"{args.prefix}.{cls}.window_tests.tsv", by_window[cls])
        write_shuffled_values(f"{args.prefix}.{cls}.shuffled_values.tsv", shuffled_rows[cls])
        stats = {"windows": window_stats[cls], "n_sites": len(rows), "n_significant": len(hits)}
        if annotation.has_validated:
            stats["n_significant_validated"] = sum(r[6] == "yes" for r in hits)
        class_stats[cls] = stats

    if not_in_table:
        examples = ", ".join(sorted(not_in_table)[:3])
        print(
            f"degrseq cleavage-stats: warning: {len(not_in_table)} sRNA(s) with real sites are "
            f"not in --srna_table and were skipped (e.g. {examples})",
            file=sys.stderr,
        )
    for name in sorted(empty_null):
        print(
            f"degrseq cleavage-stats: warning: {name} has real sites but no shuffled sites; p = 1",
            file=sys.stderr,
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "nanopare_compat": args.nanopare_compat,
        "fc_min": args.fc_min,
        "padj_max": args.padj_max,
        "rpm_min": args.rpm_min,
        "collapse_nt": args.collapse_nt,
        "classes": class_stats,
        "n_srnas_unclassified": len(unclassified),
        "n_srnas_not_in_table": len(not_in_table),
        "srnas_empty_null": sorted(empty_null),
    }


def main(argv: Sequence[str]) -> None:
    args = parse_args(argv)
    try:
        stats = run(args)
        stats_path = Path(f"{args.prefix}.cleavage-stats.stats.json")
        stats_path.write_text(json.dumps(stats, indent=2) + "\n")
    except (OSError, CleavageStatsError) as error:
        sys.exit(f"degrseq cleavage-stats: error: {error}")
