# `cap-class`: behavioural spec

Status: **approved** (18/11/2026; implemented 18/11, DL-087; decisions DL-074–DL-076, DL-080, DL-081, DL-083, DL-085, DL-086) · Replaces (nanoPARE 6264d6f): replicate sums `bedgraph_combine.py` `endClass.sh` l.180–194, trim + peak `bed_find_peaks.py --trim 1 -V pass` l.214–220, "+1 buffer" l.308–309 (compat only), uG split `bed_uug_filter.py` l.340–348
Pipeline module: `modules/local/cap_class/` (plus the bedtools/awk modules in "Pipeline context") · Audit: `nanopare_nfcore_design/answers/17_11_2026_cap-class_audit.md`

## Purpose
EndClass decides which 5′-end features of a sample group are **capped**, i.e. transcription start sites, and which are **noncapped**: decay or cleavage products. A capped RNA's 5′ end carries an extra, untemplated G added opposite the cap during reverse transcription. A feature is called capped when at least `--min_uug` of its 5P reads carry that G (uG).

`cap-class` does all numeric work of the stage for one group of replicates:
1. sums the replicates' 5P and uG tracks per nt;
2. trims each replicated region to the nt that have signal and finds its peak;
3. computes the uG fraction and splits the features into capped and noncapped.

Class labels (P, D, U, …) come afterwards from a bedtools module and do not affect the split (DL-081).

Consumers:
- `capped.bed` (columns 1–3, 6): EndMask's cap mask.
- the merged 5P sums: EndMask's dominant isoforms.
- `peaks.bed`: the class-assignment module.
- Everything else is informational.

## Inputs
| Argument | Format | Required | Notes |
|---|---|---|---|
| `--ends_plus`, `--ends_minus` | bedgraphs, one per replicate (space-separated lists, same replicate order in all four lists) | yes | 5P ends after the TSO mask (`end-counts` + mask); 0-based half-open; a row spanning n nt gives its value at each of its n positions; overlapping rows add |
| `--ug_plus`, `--ug_minus` | bedgraphs, one per replicate | yes | cap track of the same replicates. Default mode expects the default `end-counts` track (cap at the 5′-end coordinate, DL-008); compat expects the original's (G one nt upstream). The pipeline sets both from one `nanopare_compat` parameter (DL-074) |
| `--features` | BED6: chrom, start, end, name, replicate count, strand | yes | replicable regions from the `REPLICABLE_REGIONS` module |
| `--lengths` | chromosome lengths, columns 1–2 = name, length (`.fai` works) | yes | chromosome set and output order |
| `--prefix` | text | yes | output file prefix (the group name, e.g. `TG`) |

## Outputs
| File | Format | Consumer |
|---|---|---|
| `<prefix>.plus.bedgraph`, `<prefix>.minus.bedgraph` | merged 5P sums (rule 2) | EndMask (`dominant-isoforms`) |
| `<prefix>.uG.plus.bedgraph`, `<prefix>.uG.minus.bedgraph` | merged uG sums (rule 2) | none (inspection; compat golden test); always written (DL-086) |
| `<prefix>.peaks.bed` | BED6: chrom, trimmed start, trimmed end, name, peak offset (rule 4), strand; every feature that survived the trim | `CLASSIFY_FEATURES` |
| `<prefix>.capped.bed`, `<prefix>.noncapped.bed` | 10 columns (below) | `capped.bed`: EndMask cap mask (columns 1–3, 6); both: `CLASSIFY_FEATURES` |
| `<prefix>.cap-class.stats.json` | see Stats file | MultiQC |

`capped.bed` / `noncapped.bed` columns:
1. chrom;
2. start;
3. end (exclusive; compat + 1, rule 5);
4. name `5P.capped.<n>` / `5P.noncapped.<n>`, n from 1 in output order;
5. peak offset from start (compat + 1);
6. strand;
7. closest gene ID: `.` here, filled by `CLASSIFY_FEATURES`;
8. distance to it: `.` here, filled by `CLASSIFY_FEATURES`;
9. class: `.` here, filled by `CLASSIFY_FEATURES`;
10. uG fraction (rule 8).

## Rules
1. **Replicate sums.** For each strand, the per-nt sum over the replicates, separately for 5P (`--ends_*`) and uG (`--ug_*`). Default: exact sums, used unrounded by rules 3–8. Compat: each nt's sum is rounded to 2 decimals (Python `round`), and the rounded values are what rules 3–8 use, as the original reads its written files.
2. **Sum bedgraphs.**
   - Positions with sum 0 are not written. Neighbouring nt with equal values form one row.
   - Default: chromosomes in `--lengths` order, values with `--decimals` fixed decimals (DL-085).
   - Compat: chromosomes in sorted name order, values as Python repr of the 2-decimal value (`1.0`, `0.67`). A chromosome-strand whose only signal is at position 0 gets no row (quirk 37, known issue 37).
3. **Trim.** For each region in `--features`, on its own strand: the first and last nt with 5P sum > 0. Trimmed span = [first, last + 1).
   - A region with no 5P signal is dropped and counted (`dropped_zero_signal`).
4. **Peak.** The nt with the highest 5P sum in the trimmed span. Ties go to the rightmost nt, on both strands. Offset = peak − trimmed start.
   - `peaks.bed` column 5 = this offset in both modes, because `CLASSIFY_FEATURES` builds its 1-nt peaks from it.
5. **Window and written coordinates.**
   - Default: the window for rules 6–8 is the trimmed span `[start, end)`. Columns 2, 3 and 5 are written as computed. No buffer is needed: the default uG track records each cap G at the read's 5′-end coordinate, which lies inside the span (DL-074).
   - Compat: the window and written end are `end + 1`, and column 5 is offset + 1 (the original's "+1 buffer", known issue 15).
6. **Totals.** all = Σ 5P sum over the window on the feature's strand. ug = Σ uG sum over the same window and strand.
7. **Minimum.** A feature with all < `--min_reads` (1) is written to neither output and counted (`dropped_below_min_reads`; DL-075).
8. **Split.** fraction = ug / all, compared unrounded. fraction ≥ `--min_uug` → capped, otherwise noncapped (exactly `--min_uug` is capped).
   - With `--min_uug 0`, every remaining feature is capped (DL-080).
   - Compat: `--min_uug 0` selects the original's automatic cutoff instead:
     1. histogram of all fractions in 101 bins of width 0.01 (bin = floor(100 × fraction));
     2. the bins with the lowest count among bins 0–24 (fractions 0–0.25);
     3. cutoff = mean of those bins' lower edges, if several are tied.
   - Column 10 = fraction rounded to 3 decimals. Default: fixed 3 decimals (`0.545`, `0.000`). Compat: Python repr of `round(fraction, 3)` (`0.545`, `0.0`).
9. **Feature identity.** A feature is keyed by (chrom, start, end, strand). Two features with the same start but different ends or strands are distinct (known issues 52, 53).
10. **Order and names.** Rows are ordered by chromosome in `--lengths` order, then start, end, strand (`+` first). Names are numbered in this order, separately for capped and noncapped. `CLASSIFY_FEATURES` re-sorts and renames the labelled files as the original does (l.350–352).
11. Output does not depend on input row order, replicate order, `--cpus` or `PYTHONHASHSEED`.

## Pipeline context (not part of degrseq; local modules, bedtools + awk)
- **`REPLICABLE_REGIONS`** (before `cap-class`; original l.164–176, l.202–205):
  - Pool every replicate's `kde-peaks` features and add the sample as column 7.
  - `bedtools sort | bedtools merge -s -d 0 -c 6,7 -o distinct,count_distinct`, keep count ≥ 2, name `<group>.<n>` (DL-083).
  - Compat: `-c 6,5 -o distinct,count`, keep count > 1 (known issue 51).
- **`CLASSIFY_FEATURES`** (after `cap-class`; original l.224–323, l.350–352):
  1. 1-nt peaks from `peaks.bed`.
  2. The original's 8 `bedtools closest`/`intersect` calls against `class.single_exon_genes.bed`, `class.terminal_exons_by_gene.bed` and `class.exons_by_gene.bed`, with scores S 7, P 6, D 5, U 4 (≤ `upstream` nt, default 500), I 3, A 2, IA 1, N 0.
  3. Per feature key, the highest score, then mapped to its letter.
  4. Columns 7–9 joined onto `capped.bed` / `noncapped.bed` by key (default coordinates; compat end + 1 on both sides).
  5. `bedtools sort`, renamed `5P.capped.<n>` / `5P.noncapped.<n>`, and both concatenated and sorted into `<group>.5P_features.bed`.
  - The original's consecutive-row grouping (known issue 52) is not reproduced; the test set cannot differ on it.
- **Checked 18/11 on the test set** (`answers/scratch/16_11_2026_cap-class/reorder_check_18_11/`), splitting before classification with the original helper:
  - capped 323 and noncapped 279 identical to golden in coordinates, strand and uG fraction;
  - `5P_features.bed` byte-identical after joining the golden labels.

## Edge cases
- **No features, or all dropped:** empty `capped.bed`, `noncapped.bed`, `peaks.bed`; sums still written; stats with zeros; exit 0.
- **Feature on a chromosome missing from `--lengths`:** error naming it. The original drops it silently at the trim.
- **Bedgraph row past the chromosome end, or a negative or non-finite value:** error.
- **Lists of different lengths in `--ends_*` / `--ug_*`:** error.
- **A single replicate:** allowed. The sums equal that replicate (`REPLICABLE_REGIONS` normally prevents this).
- **uG without 5P at the same nt:** allowed. uG counts only through the window sums. A fraction > 1 is not possible from default `end-counts` tracks (one count per read); default mode stops with an error, which also catches a compat (per-G) track fed to default mode (DL-074). Compat accepts it, as the original.
- **Feature at position 0 or ending at the chromosome end:** valid. Compat end + 1 may equal the chromosome length + 1, as in the original.
- **Tied peak values:** the rightmost nt (rule 4).
- **Chromosome names and IDs** are literal strings (DL-064); test with `Su(var)3-9`.

## Quirks of the original
| # (known_issues) | Original behaviour | Decision | Expected diff vs golden (test set) |
|---|---|---|---|
| 15 | "+1 buffer" extends only the end; + strand cap G at `start − 1` never counted | fix + compat (DL-074) | ends and peak offsets − 1 on every feature; with the original uG tracks shifted to the 5′ end (default format): capped 323 → 330 (+ 152 → 159: 8 gained, 1 lost; − 171 unchanged; checked 18/11) |
| 50 | `--min_uug 0` → automatic cutoff | fix + compat (DL-080) | none at the default 0.1 |
| 51 | replicable = merged feature count | fix + compat (DL-083), in `REPLICABLE_REGIONS` | none (all 602 regions have ≥ 2 replicates) |
| 16 | features with < 1 read dropped from both outputs | reproduce, counted in stats (DL-075) | none |
| 52 | non-adjacent duplicate rows survive best-class selection | fix, no compat (`CLASSIFY_FEATURES`) | none (0 on the test set) |
| 53 | same-start features overwrite each other in the uG split | fix (rule 9), no compat | none |
| 37 | position-0-only chromosome-strand lost in the sums | fix + compat | none |
| — | replicate sums rounded to 2 decimals before use | default exact (DL-085); compat rounds | none on the test set (G2: 0 of 602 spans or peaks differ) |
| — | column 5 + 1 in the buffer step | default true offset; compat + 1 | peak offsets − 1 (no consumer) |
| — | `peaks.bed` column 7 = trim count | not written | none (no consumer) |
| — | zero-signal region dropped with a stdout warning | drop, count in stats | none |
| — | chromosome missing from the length table dropped silently | error | none on valid input |
| — | gene counts tables (l.396–412) | not ported (DL-076; Q19 to the authors) | files absent |
| — | `capmask` bedgraphs (l.358–386) | not ported here; EndMask masks with bedtools (13/10) | files absent |

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|
| `--min_uug` | 0.1 (number in [0, 1]; 0 = every feature capped) | `--uug` / `UUG=0.1` (l.105–107) → `bed_uug_filter.py --min_uug` |
| `--min_reads` | 1 | `bed_uug_filter.py --minimum` (default 1) |
| `--decimals` | 4 (DL-085) | 2 (fixed in `bedgraph_combine.py`) |
| `--lengths` | required | `length.table` |
| `--nanopare_compat` | off | — |
| `--prefix` | required | `$SAMPLE_TYPE` |
| (pipeline) `upstream` | 500 | `--upstream` / `UPSTREAM=500` (l.111), used by `CLASSIFY_FEATURES` |

## Stats file
`schema_version`, `mode`, `replicates` (number), `min_uug` (as used; the automatic value in compat), `min_reads`, and:
- `regions_in`;
- `dropped_zero_signal`;
- `dropped_below_min_reads`;
- `capped`, `noncapped`, per strand;
- `ug_fraction_histogram`: 101 bins of 0.01 (bin 100 = fraction 1).

## Tests
- **Unit:** one per rule and edge case, with handmade fixtures:
  - sums of 3 replicates incl. overlapping rows, both modes (0.333 + 0.334 → 0.67 in compat);
  - trim to first and last non-zero nt;
  - tied peaks → rightmost (toy T1: `[100,120)` with 3 @105, 5 @110, 5 @112 → span 105–113, offset 7);
  - + strand feature whose reads all start at its first nt and carry a G: capped in default with a default track; 0 % in compat with an original track;
  - fraction exactly at the cutoff;
  - `--min_uug 0`: all capped (default) vs automatic cutoff (compat; toy value 0.126 from the audit);
  - 0.9-read feature dropped and counted;
  - same start on both strands → two features;
  - position-0 quirk, both modes;
  - chromosome missing from `--lengths`;
  - shuffled inputs and replicate order → same output.
- **Golden (compat):**
  - inputs: nanoPARE test set, group `flower` (`flower_1–3` EndGraph `5P.*_mask` / `uG.*_mask` bedgraphs, about 2.2 MB, gzipped) plus `flower.rep.bed`;
  - expected outputs: the 4 merged bedgraphs and `peaks.bed` (columns 1–6) byte-identical to golden;
  - `capped.bed` / `noncapped.bed` equal to golden in columns 1–3, 6 and 10;
  - reference: `nanoPARE/results/EndClass/flower/` (rerun 16/11, 14/14 identical).
- **Pipeline (`CLASSIFY_FEATURES` + `cap-class` compat):** `flower.5P_features.bed` byte-identical (checked by hand 18/11).
- **Default vs golden:** the original uG tracks moved to the 5′-end coordinate (+ strand +1, − strand −1) and capped at the 5P value per nt stand in for default `end-counts` tracks (the original counts every G, so uncapped shifted tracks exceed 5P at 181–208 positions per file and default mode stops). Result (`tests/test_cap_class_golden.py`): same 602 features, every end and peak offset − 1; trimmed spans identical (G2); capped 327 (+ 157, − 170): 7 + strand features gained (quirk 15), 3 lost (one G belonging to a read outside the feature, `Ath_chr2:489670 +`; two whose golden fraction relied on several G's per read, `225715 −`, `722622 +`).

## Open
None.
- **G1, resolved 18/11:** default with shifted tracks gives 330 capped, not the audit's strand-aware 331. The difference is one feature, `Ath_chr2:489670–489672 +` (golden `5P.293`, uG 0.125): its only G sits at 489671 in flower_1, where no 5P read starts at 489672. The original window counts it; the default window does not, because the G belongs to a read starting outside the feature. Likely a G from a read missing from the 5P track (the uG leak of email Q1) [I]; default `end-counts` does not produce such G's.
- **G2, resolved 18/11:** unrounded sums change no trimmed span or peak on the test set (0 of 602).

## Decided
- D1, sum format → DL-085; D2, uG sums always written → DL-086.
- Decision 1, window → DL-074; 2, `--min_uug 0` → DL-080; 3, replicate counting → DL-083; 4, < 1 read → DL-075; 5, gene counts → DL-076; 6, one subcommand → DL-081.
