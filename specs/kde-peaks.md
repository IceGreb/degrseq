# `kde-peaks`: behavioural spec

Status: **approved** (16/11/2026, DL-072; D1–D4 decided 16/11, DL-065–DL-068; `unstranded` DL-069; G1 DL-071) · Replaces (nanoPARE 6264d6f): `bedgraph_combine.py` (subtraction) `endGraph.sh` l.255–263, `bedgraph_kernel_density.py` l.277–294, `bedgraph_thresh_to_bed.py` l.296–312, strand merge + `bedtools sort` l.316–323, `bed_find_peaks.py -V rpm` l.325–331, RPM filter l.333; bandwidth cap l.249–253
Pipeline module: `modules/local/kde_peaks/` · Audit: `nanopare_nfcore_design/answers/12_11_2026_kde-peaks_audit.md`

## Purpose
EndGraph finds the 5′-end features of one 5P library:
1. Subtract the BODY 5′ ends, which estimate the ends produced by ordinary fragmentation, from the 5P ends scaled by the sample's factor (`scale-factor`).
2. Smooth the difference with a Laplace kernel.
3. Call every run of positive density at least `--min_length` nt long a feature.
4. Report each feature's strongest 5P position and its RPM, and keep features with at least `--min_rpm`.

Consumer: EndClass (pipeline replicate merge, then `cap-class`). It uses only the feature coordinates and strand; the other columns are informational.

## Inputs
| Argument | Format | Required | Notes |
|---|---|---|---|
| `--ends_plus`, `--ends_minus` | bedgraph, 5P ends after the TSO mask | yes | 0-based half-open; a row spanning n nt gives its value at each of its n positions; rows may be merged or per-nt; overlapping rows add |
| `--body_plus`, `--body_minus` | bedgraph, BODY ends | yes | strand as aligned (`end-counts --library BODY`) |
| `--scale_factor` | text file, one number (`<prefix>.scale_factor.txt` from `scale-factor`) | yes | used as given (the cap is applied by `scale-factor`) |
| `--lengths` | chromosome lengths: tab-separated, columns 1–2 = name, length (`.fai` works) | yes | defines the chromosome set and the output chromosome order |

## Outputs
| File | Format | Consumer |
|---|---|---|
| `<prefix>.end_features.bed` | 9 tab-separated columns, no header (below) | EndClass replicate merge (columns 1–3, 6) |
| `<prefix>.density_plus.bedgraph`, `<prefix>.density_minus.bedgraph` | only with `--write_density`: smoothed density, positive values, rows of equal values merged | none (inspection; the original always writes them) |
| `<prefix>.kde-peaks.stats.json` | see Stats file | MultiQC |

`end_features.bed` columns:
1. chrom;
2. start;
3. end (exclusive; compat: inclusive last position, quirk 12);
4. name `5P.<plus|minus>.<chrom>.<n>`;
5. dominant 5P position, relative to start (rule 8);
6. strand `+`/`-`;
7. RPM (rule 9);
8. highest density peak, relative to start (rule 7);
9. up to 9 further density peaks, comma-separated, or `-`.

## Rules
1. **Pairing.** 5P plus is paired with BODY minus and 5P minus with BODY plus (`--body_strand reverse`, default); with `forward`, plus with plus and minus with minus; with `unstranded`, both 5P strands are paired with BODY plus + BODY minus summed per nt (DL-069). Same parameter and meaning as `scale-factor`.
2. **Difference.** Per chromosome, strand and nt: d = S × E − B, where S is the scale factor, E the 5P value and B the paired BODY value (0 where absent). Default: exact (no rounding). Compat: each d is rounded to 2 decimals (Python `round`) before rule 4, and a chromosome-strand whose only non-zero d is at position 0 contributes nothing (quirk 37).
3. **Kernel.** W = `--kernel_extent` × H, where H = `--bandwidth`. Weights k(x) = exp(−|x| / H) for integer x = −W … +W, normalised to sum 1 over these 2W + 1 values. Default: all 2W + 1 weights are used. Compat: the weight at x = +W is not used (the others are not renormalised; quirk 10). Precision: see G1.
4. **Density.** D(p) = Σₓ d(x) · k(p − x) over all nt x of the same chromosome and strand. Negative d lowers the density. D is computed for 0 ≤ p < chromosome length. Compat: for p ≥ 0 without the upper bound (quirk 11). Positions p < 0 are never reported (both modes).
5. **Rounding of the density.** D is rounded to `--decimals` (3) decimals. A position is **positive** if its rounded D > 0. The rounded values are what rules 6–7 use and what `--write_density` writes. This rounding is kept in default mode because it sets the detection limit: without it, every feature would extend to the kernel's full reach (DL-066).
6. **Features.** A feature is a maximal run of consecutive positive positions on one chromosome and strand with length ≥ `--min_length`. Start = its first position. End = last + 1 (compat: last; quirk 12). Names: n counts this strand's features on this chromosome from 1 in ascending start order, before RPM filtering, so kept features may have gaps in n.
7. **Density peaks (columns 8–9).**
   - Walk the feature's rounded densities from its first position, with the value before the first position taken as 0.
   - A peak is the last position of a run of equal values that was entered by a rise and is followed by a lower value.
   - If the feature has no such peak and its final run was entered by a rise, the final position is the only peak.
   - Peaks are ordered by density descending, ties by position descending. Column 8 = the first; column 9 = the next ≤ 9 joined by `,`, or `-` if there are none.
8. **Dominant 5P position (column 5).** Over the feature's positions [start, end), the position with the largest E (5P value on the feature's strand). Ties go to the largest position; if all E = 0, the last position. Compat: over [start, last), i.e. the feature's last nt is never considered (follows from quirk 12).
9. **RPM (column 7).** RPM = (Σ E over the positions of rule 8) / N × 10⁶. Default: N = total of all 5P values on both strands (`math.fsum`), and RPM is written with 2 fixed decimals. Compat: N = int(2 × total of `--ends_plus`) (quirk 13), and RPM is written as Python's repr of `round(RPM, 2)` (`5.14`, `0.0`).
10. **Filter.** Keep features with RPM ≥ `--min_rpm`. Default compares the exact RPM; compat compares the rounded value (as the original's awk on column 7).
11. **Order.** Chromosomes in `--lengths` order; within a chromosome, by start, then end, then strand (`+` first). Compat: same, except that rows with equal chromosome and start keep no defined order (G2).
12. **Number format.** Default: column 7 with 2 fixed decimals; density values with `--decimals` fixed decimals. Compat: Python repr of the rounded value (as the original). Integers (columns 2, 3, 5, 8, 9) are plain.
13. Output does not depend on input row order, `--cpus` or `PYTHONHASHSEED`.

## Edge cases
- **Chromosome in a bedgraph but not in `--lengths`:** error naming the chromosome, in both modes. The original processes it in the smoothing step, crashes in the threshold step, and drops its features silently in the peak step.
- **Bedgraph row past the chromosome end, negative or non-finite value:** error.
- **Scale factor missing, ≤ 0 or non-finite:** error.
- **No 5P values (N = 0):** there can be no positive density, so no features. Write an empty BED and stats; not an error. Compat with an empty `--ends_plus` but features on minus: error, because the original divides by zero.
- **No BODY values:** allowed; d = S × E.
- **Features within W nt of position 0:** density before 0 is lost (both modes); the feature starts at 0.
- **Density near the chromosome end:** default stops at length − 1; compat writes past it.
- **A feature of exactly `--min_length` nt is kept; one nt shorter is dropped.**
- **Equal density peaks or equal 5P maxima:** the larger position wins (rules 7, 8).
- **Chromosome names** are used as literal strings (DL-064). Test with a name containing regex metacharacters, e.g. `Su(var)3-9`.

## Quirks of the original
| # (known_issues) | Original behaviour | Decision | Expected diff vs golden (test set) |
|---|---|---|---|
| 10 | kernel weight +W dropped | fix + compat (DL-065) | 667 / 1,093 flower_1 features shift (mostly end + 1); 6–9 fewer features per sample; replicable 602 → 598 |
| 11 | density past the chromosome end; `-L` ignored | fix + compat (DL-065) | none (no feature near the end) |
| 12 | BED end = last position; last nt left out of RPM and dominant position | fix + compat (DL-065) | all features + 1 nt; kept + 5 / + 8 / + 5; replicable 602 → 607 |
| 13 | RPM denominator = int(2 × plus total) | fix + compat (DL-065) | RPM values change; kept set unchanged (the filter removes only 0-read features here) |
| 37 | subtraction: a chromosome-strand with its only value at position 0 gets no row | fix + compat | none |
| — | difference rounded to 2 decimals before smoothing | default exact; compat rounds (DL-066) | G3 |
| — | chromosome missing from the length table: crash or silent drop | error in both modes | none on valid input |
| 14 | Gaussian kernel width = H√2 | not carried over: Laplace only (DL-068) | none (Laplace is the default) |
| 47 | BODY strand flip hard-coded | `--body_strand` `reverse` (default, = original) / `forward` / `unstranded` | none |
| — | bandwidth > 30 silently capped to 30 | error for bandwidth > 30 | none at defaults |
| — | README says default RPM 5 "for test data"; script default 0.5 | 0.5 (script, golden run) | none |
| — | number repr (`14.0`, `5.14`) | fixed decimals in default; repr in compat | text only |

With all of quirks 10–13 fixed together: replicable features 602 → 606 (+0.7 %).

## Parameters
| Flag | Default | Original equivalent |
|---|---|---|
| `--bandwidth` | 15 (integer 1–30) | `--bandwidth` / `BANDWIDTH`; cap `BANDWIDTH_CAP=30` |
| `--kernel_extent` | 3 (integer 1–6) | `bedgraph_kernel_density.py -S 3` |
| `--decimals` | 3 | `-D 3` |
| `--min_length` | 10 | `bedgraph_thresh_to_bed.py -M 10` |
| `--min_rpm` | 0.5 | `--rpm` / `RPM=0.5` |
| `--body_strand` | `reverse` (`reverse` / `forward` / `unstranded`) | `flip_5p="true"` (hard-coded) |
| `--write_density` | off | (always written: `5P.<strand>_smooth.bedgraph`) |
| `--nanopare_compat` | off | — |
| `--prefix` | required | `$SAMPLE_NAME` |
| (not exposed) threshold | > 0 | `-T 0` |
| (not exposed) kernel | Laplace | `--kernel laplace` (Gaussian dropped, DL-068) |

## Stats file
`schema_version`, `mode` (`default`/`nanopare_compat`), `scale_factor`, `body_strand`, `bandwidth`, `kernel_extent`, `decimals`, `min_length`, `min_rpm`, `total_5p` (N as used), and per strand (`plus`, `minus`):
- `positive_positions`;
- `features_called` (rule 6);
- `features_kept` (rule 10);
- `features_without_reads` (Σ E = 0).

## Tests
- **Unit:** one per rule and edge case, with handmade fixtures:
  - unit point → kernel shape, sum and symmetry; compat asymmetry, with the toy values 0.174 / 0.186 at ±44–45 for v = 100;
  - negative difference lowers density;
  - point near 0 and near the chromosome end, both modes;
  - 9 vs 10 nt features;
  - plateau and tie peaks (toy vectors `[1,1,3,1,1,3,1,1,1,1]` → 5 then `2`; `[2,2,2,5,5,5,1,1,1,1]` → 5, `-`; `[1,3,3,1,3,3]` → 2);
  - all-zero and tied 5P maxima;
  - RPM denominator both modes;
  - `--body_strand` all three values (unstranded: BODY on either strand lowers both strands' density);
  - chromosome missing from lengths;
  - input row order shuffled → same output.
- **Golden (compat), byte-identical:** flower_1 of the nanoPARE test set, with inputs `5P.{plus,minus}_mask.bedgraph`, `BODY.{plus,minus}.bedgraph`, factor 0.3924, `resources/length.table`. Expected outputs: `flower_1.end_features.bed`, and with `--write_density` the two `5P.*_smooth.bedgraph`. Fixture ≈ 4.3 MB.
  - The ordering exception for rows with equal chromosome and start (G2): there is 1 such pair in flower_1's golden output.
  - flower_2 / flower_3 (factors 0.2100 / 0.2121) are checked in the audit scratch run, not as fixtures (size).
- **Default vs golden:** the difference must be exactly the quirk table: the 5 kept-feature gains from quirk 12, the shifts from quirk 10, and G3's measured effect.

## Open
- **G1, kernel precision (procedure decided, DL-071).** The original's weights are single precision; double precision differs by ≤ 2.7 × 10⁻⁹ per weight, which can flip a density's 3rd decimal only next to a rounding boundary. Implement in double precision; if the compat golden test fails only through such flips, compat uses single-precision weights and this item records it. Does not block approval.
- **G2, tie order.** The original orders rows with equal chromosome and start by `bedtools sort`, which leaves their order undefined: the 3 pre-filter ties on the test set come out in 3 different patterns. Proposal: compat uses rule 11's order, and the golden test compares after sorting both files by chromosome, start and strand (2 golden rows affected).
- **G3, effect of dropping the 2-decimal rounding** (rule 2 default). Not measured; to be measured when implementing (features called, kept and replicable on the test set) and recorded in the quirk table.

## Decided
- `--body_strand unstranded` → DL-069; G1 procedure → DL-071.
- D1 → DL-065 (fix quirks 10–13, 37 with compat); D2 → DL-066 (no intermediate rounding in default; density rounding kept as the detection limit); D3 → DL-067 (9 columns, density bedgraphs only with `--write_density`); D4 → DL-068 (Laplace only).
