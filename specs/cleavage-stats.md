# `cleavage-stats`: behavioural spec

Status: **approved 30/10/2026** (DL-052; drafted 13/10, D1–D4 and G3 decided 16/10, D5 30/10, amended from the 30/10 audit of the output code) · **implemented 30/10/2026** · Replaces (nanoPARE 6264d6f): `EndCut.pipeline/scripts/EndCut_step2.v03.git.R` (385 lines; R packages `survcomp`, `stringr`)
Pipeline module: `modules/local/cleavage_stats/` · Design-note section: §6 EndCut

## Purpose
Decides which predicted sRNA target sites show significant cleavage. For each real sRNA it
compares the fold change (5′-end signal at the site vs its flanks) and the Allen score
(pairing quality) of its detected sites with the same quantities for that sRNA's 1,000
shuffled versions. The two empirical p-values are combined (Fisher), BH-adjusted across sites,
and filtered. This is the step that outputs cleavage sites, the pipeline's end goal.

## Inputs
| Argument | Format | Required | Notes |
|---|---|---|---|
| `--real_sites_w20`, `--real_sites_w50` | EndCut step-1 table for the real sRNA set, WIN = 20 and 50 | yes | TSV with header; columns used by name: sRNA, Transcript, Slice.site (0-based), AllenScore, flankMax.rpm, slice.site.rpm. Fields are whitespace-stripped. Any table may be gzip-compressed (`.gz`) |
| `--shuffled_sites_w20`, `--shuffled_sites_w50` | the same for all shuffled sets, pooled in one table, no set column (DL-050); rows that repeat across sets are separate observations and must not be de-duplicated | yes | Original: 1,000 files per window, `anno.mir.tas.<i>.shuffled.fa…_detectedSites.txt` |
| `--srna_table` | TSV: `sRNA` (required), `class` (`miRNA` / `tasiRNA`, optional), `family` (optional) | yes | Original: `Ath_annotations/anno.{miRNA,tasiRNA}.tsv`, class and family from the name only. Missing class/family → name rule (Q1); unclassified sRNAs skipped with a warning (D1, decided) |
| `--precursors` | TSV with header: `sRNA`, `precursor` (gene ID) | optional | Original: `miRBase21_and_TAIR10_miRNA_tasiRNA_precursors.txt`; used to drop self:self hits |
| `--validated_targets` | TSV with header: `sRNA`, `targetAGI` (validated target gene, 5′ RACE) | optional | Original: `sRNA_targets_5pRACE.tsv` (same columns); only annotates |
| `--gene_aliases` | TSV with header: `locus_name` (gene ID), `symbol`; the first row per gene is used | optional | Original: `gene_aliases_20130130.txt` (TAIR, same columns); only annotates |

## Outputs
| File | Format | Consumer |
|---|---|---|
| `<prefix>.<class>.detected_sites.tsv` | One row per site; header as the original writes it: target.transcript, sRNA, slice.site, Allen.score, slice.site.rpm, fold.change, verified. (R renames `verified?`), common.name, interaction, adj.p.val. Sorted by rule 12; numbers as rule 14. `<class>` = `miRNA` / `tasiRNA`; always written, header only if no sites | user / report |
| `<prefix>.<class>.significant_sites.tsv` | Rows of `detected_sites.tsv` passing rule 9 on the reported values; same columns and order. New (the original computes but doesn't write them) | user / report |
| `<prefix>.cleavage-stats.stats.json` | see Stats | MultiQC |

## Rules (per sRNA class, per window W ∈ {20, 50})
1. **Fold change:** `fc = (slice.site.rpm + 1) / (flankMax.rpm + 1)`, for real and shuffled rows.
2. **Real sites, filter:** keep rows whose sRNA is in the class table. sRNAs with real sites that aren't in `--srna_table` are skipped with one warning and counted (`n_srnas_not_in_table`), so a naming mismatch is visible.
3. **Self:self removal (real only):** drop a row if the transcript's gene ID (transcript ID up
   to the first `.`) equals the sRNA's precursor gene. An sRNA without a precursor keeps its rows (the original stops with an error). An sRNA may have several precursor rows (e.g. miRNAs from several loci); a row is dropped if its gene matches any of them (the original stops on duplicate rows).
4. **Redundant-site collapse (real only):** sort rows by Allen score ascending, then
   slice.site.rpm descending, then sRNA ascending (intended order; see G1). Per transcript, keep the first
   row, then keep each later row only if its slice site is more than 5 nt from every row
   already kept for that transcript.
5. **Shuffled sites:** all sets pooled per sRNA. No self:self or redundancy filtering.
   Missing values (empty, `NA`, `NaN`) become 0, and fc is 0 if rpm or flankMax.rpm is missing; the slice site isn't used. Missing values in real rows are an error. All values must be finite and ≥ 0.
6. **Empirical p-values, per real row, against that sRNA's pooled shuffled rows:**
   `p_fc = 1 − F_fc(fc)` and `p_allen = F_allen(Allen)`, where `F(q)` = fraction of shuffled values ≤ q (R `ecdf`).
   If the sRNA has no shuffled rows, `p_fc = p_allen = 1` and a warning is written (D3).
7. **Fisher combination:** `X = −2(ln p_fc + ln p_allen)`, `p = P(χ²₄ ≥ X)`. Any p = 0 gives combined p = 0; silence numpy's divide-by-zero warning (G2).
8. **BH adjustment** over all real rows of the class and window (`p.adjust(method='fdr')`).
9. **Significant:** `fc > 1.0` and `adj_p < 0.05` and `slice.site.rpm > 0.1`.
10. **Combining windows:** outer-join W20 and W50 results on `Transcript_Slice.site`. Per site,
    report the window with the smaller adjusted p-value, along with that window's sRNA, Allen,
    rpm, fc and slice site. Ties → W20.
11. **Labels:** `common.name` from the aliases (else `NA`). `verified?` = `yes` if a validated
    target with the same gene exists whose sRNA name, or family, equals the sRNA's family; a validated sRNA's family is its `--srna_table` family if listed, else the class's name rule (rule 13). Matching the name too lets `miR2934` validate `miR2934-3p` with family `miR2934`; with the name rule alone (the original) it can't change a result. `interaction` = `sRNA:common.name`, falling back to `sRNA:Transcript`, with the reported sRNA.
    When joining windows (rule 10), each label takes the non-missing value; `verified?` = `yes` if either window says yes (D2). Compat: R's `min()` (Q3): `verified?` = `no` if either window says no; `common.name` = `NA` for a site seen in one window only; `interaction` uses the W20 sRNA when W20 has the site (Q9).
12. **Output order:** adj p ↑, fc ↓, sRNA ↑, Transcript ↑, Slice.site ↑ (G3). Compat: sites ordered by the text `Transcript_Slice.site` (byte order; R `merge` sorts on it), then a stable sort on adj p.
13. **Class and family (D1):** from `--srna_table` where given; otherwise `miR…` → miRNA, family = first 6 characters; `TAS…` → tasiRNA, family = first 4. Anything else is unclassified: skipped, warned, counted. A class other than `miRNA` / `tasiRNA`, or a duplicate sRNA, is an error.
14. **Numbers** are written as R's `write.table` does: up to 15 significant digits, and per value scientific notation (`1e-04`) when it is shorter than fixed notation. Integers (slice site) as is.

## Edge cases
- An sRNA with real sites but no shuffled sites → `F` undefined (R `ecdf` errors on empty input, aborting the whole script). **Decided (D3):** p = 1, warning, listed in stats `srnas_empty_null`.
- A sRNA class with no detected sites → empty outputs, no error.
- Ties in `fc` or Allen at the ECDF boundary are handled by rule 6's "≤".
- A site detected in only one window → its values come from that window (rule 10 outer join). The original stops with an error when a site is only in W50 or a window has no sites (Q10).

## Quirks of the original
| # | Original behaviour | Decision | Expected diff |
|---|---|---|---|
| Q1 | Class and family come from name prefixes and fixed character counts (`miR…` → 6 chars, `TAS…` → 4) | **fix (D1)**: explicit `--srna_table` columns; prefix rule as fallback (rule 13) | none on Arabidopsis without class/family columns |
| Q2 | Arabidopsis-specific inputs and paths are hard-coded (`Ath_annotations/…`, TAIR aliases) | **fix**: all optional arguments | — |
| Q3 | `print.tables` takes `min()` of the two `verified?` strings and the two `common.name` strings (lexicographic), and calls `min(…, na.rm=TRUE)` with the `na.rm` misplaced on `common.name` | **fix + compat (D2)**: default = non-missing value ("yes" if either window says yes); compat reproduces `min()`. [V] 15/10: `min("ARF17", NA)` = `NA` in R 4.5.3, so a site detected in only one window loses its `common.name` (l.336) | labels only; calls unchanged |
| Q4 | The two windows are run with a hard-coded 1,000 shuffles (`get.tars`), and `shuff.num` is not taken from the data | **moot (DL-050)**: one pooled table, no set number; any number of sets works | none |
| Q5 | Per-sRNA ECDF PDFs are always written | **dropped from the package (D4, 30/10)**: plots are left to the pipeline, drawn from the output tables (the per-sRNA ECDFs would also need the shuffled values, which are not written) | no PDFs |
| Q7 | An sRNA with an empty null aborts the whole script | **fix (D3)**: p = 1 + warning | runs that crashed now complete |
| Q8 | Output sorted by adj p only (stable on input order) | **fix + compat (G3)**: full tie-break (rule 12) | row order among ties |
| Q6 | The p = 0 Fisher case (rule 7) | must match | G2 |
| Q9 | `interaction` takes the W20 sRNA whenever W20 has the site, even when the W50 row is reported | **fix + compat**: default uses the reported sRNA | `interaction` only, when the two windows' collapses keep different sRNAs |
| Q10 | `print.tables` stops with an error (`if (NA == …)`) for a site seen only in W50, and `get.tars` builds a table without `.x`/`.y` columns when a window has no sites, which also stops it [V 30/10, R 4.5.3] | **fix, both modes**: outer join (rule 10) | runs that stopped now complete; no golden output to match |

## Parameters
| Flag | Default | Original |
|---|---|---|
| `--fc_min` | 1.0 | hard-coded |
| `--padj_max` | 0.05 | hard-coded |
| `--rpm_min` | 0.1 | hard-coded |
| `--collapse_nt` | 5 (≥ 0) | hard-coded |
| `--nanopare_compat` | off | — |
| `--prefix` | required | — |

## Stats file
`schema_version`, the thresholds and `nanopare_compat`; per class: per window `n_real_rows` (in the class), `n_after_self_self`, `n_after_collapse`, `n_srnas_tested`, `n_shuffled_rows`; then `n_sites`, `n_significant`, and `n_significant_validated` (only with `--validated_targets`). Overall: `n_srnas_unclassified` (D1), `n_srnas_not_in_table` (rule 2), `srnas_empty_null` (sorted names, D3). Shuffled sets can't be counted (DL-050).

## Tests
- Unit: each rule on handmade tables (fc pseudocount; ±5-nt collapse ordering; ECDF "≤"; Fisher incl. p = 0; BH vs `p.adjust` reference values; window join chooses min adj p; label fallbacks; D1 table vs name rule incl. `ame-miR-1` unclassified and `miR2934-3p` with an explicit family; D2 one-window site keeps its name (compat: `NA`); D3 empty null → p = 1 + warning; G3 tie order).
- Golden: the golden run (Phase 5b, 29/10) of `EndCut_step2.v03.git.R` + P3 on the test set: flower_1–3, both classes, with all 1,000 shuffled step-1 tables pooled per window. Compat output must be byte-identical (`tests/fixtures/cleavage_stats_golden/`). A scratch re-implementation already matched all 41 adjusted p-values (`package_build/30_10_2026_cleavage-stats_pvalue_test_run.md`).

## Open
- **D5, decided 30/10 (DL-051, provisional):** keep the original rules 6–7 (no +1, fold-change ties not counted, Fisher 0). Test run: the corrected formula `(b + 1)/(n + 1)` calls 0 of the 9 golden sites, miR173 → TAS1C/TAS2 included, because nulls have 3–16 rows. Authors asked (Q14). Evidence: `nanopare_nfcore_design/package_build/30_10_2026_cleavage-stats_zero_p_values.md`, `30_10_2026_cleavage-stats_pvalue_test_run.md`.
- **G1, resolved 15/10/2026 [V]:** in R 4.5.3, `order(…, decreasing=c(FALSE,TRUE,FALSE))` (R l.79) **errors** (`'decreasing' of length > 1 is only for method = "radix"`), because `sRNA` is character (`stringsAsFactors=FALSE`). With `method="radix"` or a factor key it gives the intended order: Allen ↑, rpm ↓, sRNA ↑. Rule 4 implements that order. nanoPARE pins R 3.4.0 (not tested); the shipped outputs must have come from an environment or script version where this line ran.
- **G2, resolved 15/10/2026 [V]:** survcomp's Fisher method is `pchisq(-2*sum(log(p)), df=2k, lower.tail=FALSE)` (bhklab/survcomp `R/combine.test.R`; weights ignored). `scipy.stats.combine_pvalues(method='fisher')` (scipy 1.18.1) equals base-R values to < 1e-12 relative on 8 cases, including p = 0 → 0, (1, 0) → 0 and 1e-300 underflow → 0. `false_discovery_control(method='bh')` equals `p.adjust('fdr')`.
- **G3, decided 16/10/2026:** many sites tie at adj p = 0; rule 12 gives the tie-break. R's `order` is stable, so in compat mode ties keep input order.
- **D1–D4, decided 16/10/2026** (decision log DL-018–DL-021): D1 table + name-rule fallback (rule 13); D2 non-missing label, compat `min()` (rule 11); D3 p = 1 + warning (rule 6); D4 `--plots` optional extra (superseded 30/10, DL-054: no plots in the package).
- **adj p precision (30/10, DL-052):** compat matches the original's adj p to the 15th significant digit at worst (5 of 41 golden values differ in that digit; scipy vs R chi-square tail). Not a concern; the golden test uses a relative tolerance of 1e-12.
- Golden p-values: resolved by the golden run (29/10), which regenerated the 1,000 shuffled tables. The authors' own tables (Q3, DL-023) are no longer needed.
- **Both windows hold the same sites (checked 30/10):** EndCut step 1 builds the W20 and W50 tables from the same predicted sites; only flankMax.rpm (so fc) differs. On the golden run the real tables of all 3 samples have identical rows apart from that column. So with pipeline input, rule 4 keeps the same sRNA per site in both windows, and Q9, Q10 and the one-window or disagreeing-label cases of rule 11 (D2/Q3) can't occur; they guard against tables from different runs. What does matter: which window is reported (rule 10, ties → W20) and the output order (rule 12).
- **D4 revised 30/10 (DL-054):** the package writes tables only; any plots belong to the pipeline. No matplotlib extra.
