# `cleavage-stats`: behavioural spec

Status: **draft** (13/10/2026; D1–D4 and G3 decided 16/10/2026, see `package_build/degrseq_decision_log.md` DL-018–DL-022; awaiting approval of the whole spec) · Replaces (nanoPARE 6264d6f): `EndCut.pipeline/scripts/EndCut_step2.v03.git.R` (385 lines; R packages `survcomp`, `stringr`)
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
| `--real_sites_w20`, `--real_sites_w50` | `cleavage-window` output for the real sRNA set, WIN = 20 and 50 | yes | Columns: sRNA, Transcript, Slice.site (0-based), AllenScore, flankMax.rpm, slice.site.rpm, slice.site.hit.norm |
| `--shuffled_sites_w20`, `--shuffled_sites_w50` | the same for all shuffled sets, pooled in one table, no set column (DL-050); rows that repeat across sets are separate observations and must not be de-duplicated | yes | Original: 1,000 files per window, `anno.mir.tas.<i>.shuffled.fa…_detectedSites.txt` |
| `--srna_table` | TSV: `sRNA` (required), `class` (`miRNA` / `tasiRNA`, optional), `family` (optional) | yes | Original: `Ath_annotations/anno.{miRNA,tasiRNA}.tsv`, class and family from the name only. Missing class/family → name rule (Q1); unclassified sRNAs skipped with a warning (D1, decided) |
| `--precursors` | sRNA → precursor gene ID | optional | Original: `miRBase21_and_TAIR10_miRNA_tasiRNA_precursors.txt`; used to drop self:self hits |
| `--validated_targets` | sRNA → validated target gene (5′ RACE) | optional | Original: `sRNA_targets_5pRACE.tsv`; only annotates |
| `--gene_aliases` | gene ID → common name | optional | Original: `gene_aliases_20130130.txt` (TAIR); only annotates |

## Outputs
| File | Format | Consumer |
|---|---|---|
| `<prefix>.<class>.detected_sites.tsv` | One row per site; columns as the original `print.tables`: target.transcript, sRNA, slice.site, Allen.score, slice.site.rpm, fold.change, verified?, common.name, interaction, adj.p.val. Sorted by rule 12 | user / report |
| `<prefix>.<class>.<sRNA>.ecdf.pdf` | Per-sRNA ECDF plots; only with `--plots` (D4) | user |
| `<prefix>.<class>.significant_sites.tsv` | Rows passing the filter (rule 9) | user / report |
| `<prefix>.cleavage-stats.stats.json` | see Stats | MultiQC |

## Rules (per sRNA class, per window W ∈ {20, 50})
1. **Fold change:** `fc = (slice.site.rpm + 1) / (flankMax.rpm + 1)`, for real and shuffled rows.
2. **Real sites, filter:** keep rows whose sRNA is in the class table.
3. **Self:self removal (real only):** drop a row if the transcript's gene ID (transcript ID up
   to the first `.`) equals the sRNA's precursor gene.
4. **Redundant-site collapse (real only):** sort rows by Allen score ascending, then
   slice.site.rpm descending, then sRNA ascending (intended order; see G1). Per transcript, keep the first
   row, then keep each later row only if its slice site is more than 5 nt from every row
   already kept for that transcript.
5. **Shuffled sites:** all sets pooled per sRNA. No self:self or redundancy filtering.
   Missing values become 0.
6. **Empirical p-values, per real row, against that sRNA's pooled shuffled rows:**
   `p_fc = 1 − F_fc(fc)` and `p_allen = F_allen(Allen)`, where `F(q)` = fraction of shuffled values ≤ q (R `ecdf`).
   If the sRNA has no shuffled rows, `p_fc = p_allen = 1` and a warning is written (D3).
7. **Fisher combination:** `X = −2(ln p_fc + ln p_allen)`, `p = P(χ²₄ ≥ X)`. Any p = 0 gives combined p = 0; silence numpy's divide-by-zero warning (G2).
8. **BH adjustment** over all real rows of the class and window (`p.adjust(method='fdr')`).
9. **Significant:** `fc > 1.0` and `adj_p < 0.05` and `slice.site.rpm > 0.1`.
10. **Combining windows:** outer-join W20 and W50 results on `Transcript_Slice.site`. Per site,
    report the window with the smaller adjusted p-value, along with that window's sRNA, Allen,
    rpm, fc and slice site.
11. **Labels:** `common.name` from the aliases (else `NA`). `verified?` = `yes` if a validated
    target with the same gene and the same sRNA family exists. `interaction` = `sRNA:common.name`, falling back to `sRNA:Transcript`.
    When joining windows (rule 10), each label takes the non-missing value; `verified?` = `yes` if either window says yes (D2). Compat: R's `min()` (Q3).
12. **Output order:** adj p ↑, fc ↓, sRNA ↑, Transcript ↑, Slice.site ↑ (G3). Compat: stable sort on adj p only, input order kept.
13. **Class and family (D1):** from `--srna_table` where given; otherwise `miR…` → miRNA, family = first 6 characters; `TAS…` → tasiRNA, family = first 4. Anything else is unclassified: skipped, warned, counted.

## Edge cases
- An sRNA with real sites but no shuffled sites → `F` undefined (R `ecdf` errors on empty input, aborting the whole script). **Decided (D3):** p = 1, warning, listed in stats `srnas_empty_null`.
- A sRNA class with no detected sites → empty outputs, no error.
- Ties in `fc` or Allen at the ECDF boundary are handled by rule 6's "≤".
- A site detected in only one window → its values come from that window (rule 10 outer join).

## Quirks of the original
| # | Original behaviour | Decision | Expected diff |
|---|---|---|---|
| Q1 | Class and family come from name prefixes and fixed character counts (`miR…` → 6 chars, `TAS…` → 4) | **fix (D1)**: explicit `--srna_table` columns; prefix rule as fallback (rule 13) | none on Arabidopsis without class/family columns |
| Q2 | Arabidopsis-specific inputs and paths are hard-coded (`Ath_annotations/…`, TAIR aliases) | **fix**: all optional arguments | — |
| Q3 | `print.tables` takes `min()` of the two `verified?` strings and the two `common.name` strings (lexicographic), and calls `min(…, na.rm=TRUE)` with the `na.rm` misplaced on `common.name` | **fix + compat (D2)**: default = non-missing value ("yes" if either window says yes); compat reproduces `min()`. [V] 15/10: `min("ARF17", NA)` = `NA` in R 4.5.3, so a site detected in only one window loses its `common.name` (l.336) | labels only; calls unchanged |
| Q4 | The two windows are run with a hard-coded 1,000 shuffles (`get.tars`), and `shuff.num` is not taken from the data | **fix**: count the sets present | none |
| Q5 | Per-sRNA ECDF PDFs are always written | **fix (D4)**: optional (`--plots`, matplotlib extra `degrseq[plots]`, imported only when used) | no PDFs by default |
| Q7 | An sRNA with an empty null aborts the whole script | **fix (D3)**: p = 1 + warning | runs that crashed now complete |
| Q8 | Output sorted by adj p only (stable on input order) | **fix + compat (G3)**: full tie-break (rule 12) | row order among ties |
| Q6 | The p = 0 Fisher case (rule 7) | must match | G2 |

## Parameters
| Flag | Default | Original |
|---|---|---|
| `--fc_min` | 1.0 | hard-coded |
| `--padj_max` | 0.05 | hard-coded |
| `--rpm_min` | 0.1 | hard-coded |
| `--collapse_nt` | 5 | hard-coded |
| `--windows` | 20,50 | hard-coded |
| `--plots` | off | always on |
| `--nanopare_compat` | off | — |

## Stats file
`schema_version`; per class: n sRNAs tested, n real sites (before and after self:self and collapse), n shuffled sets found, n significant sites, n significant that are 5′-RACE validated (if `--validated_targets` given). Overall: `n_srnas_unclassified` (D1), `srnas_empty_null` (sorted names, D3).

## Tests
- Unit: each rule on handmade tables (fc pseudocount; ±5-nt collapse ordering; ECDF "≤"; Fisher incl. p = 0; BH vs `p.adjust` reference values; window join chooses min adj p; label fallbacks; D1 table vs name rule incl. `ame-miR-1` unclassified and `miR2934-3p` with an explicit family; D2 one-window site keeps its name (compat: `NA`); D3 empty null → p = 1 + warning; G3 tie order; `--plots` without matplotlib installed fails with a clear message).
- Golden: the clone ships the *final* output (`resources/detected.sites.capmasked/fb_1.anno.{mir,tas}.detected.sites.tsv`) and the real-set W20/W50 inputs, but **not** the 1,000 shuffled tables. So p-values can't be reproduced from the clone alone. Options: (a) ask the authors for the shuffled tables; (b) run the original R on a full regenerated null (needs GSTAr ×1001; Phase 5b, needs your go-ahead). Until then, test rules 1–4 and 10–11 against the shipped final table where p-values aren't involved.

## Open
- **D5, decided 30/10 (DL-051, provisional):** keep the original rules 6–7 (no +1, fold-change ties not counted, Fisher 0). Test run: the corrected formula `(b + 1)/(n + 1)` calls 0 of the 9 golden sites, miR173 → TAS1C/TAS2 included, because nulls have 3–16 rows. Authors asked (Q14). Evidence: `nanopare_nfcore_design/package_build/30_10_2026_cleavage-stats_zero_p_values.md`, `30_10_2026_cleavage-stats_pvalue_test_run.md`.
- **G1, resolved 15/10/2026 [V]:** in R 4.5.3, `order(…, decreasing=c(FALSE,TRUE,FALSE))` (R l.79) **errors** (`'decreasing' of length > 1 is only for method = "radix"`), because `sRNA` is character (`stringsAsFactors=FALSE`). With `method="radix"` or a factor key it gives the intended order: Allen ↑, rpm ↓, sRNA ↑. Rule 4 implements that order. nanoPARE pins R 3.4.0 (not tested); the shipped outputs must have come from an environment or script version where this line ran.
- **G2, resolved 15/10/2026 [V]:** survcomp's Fisher method is `pchisq(-2*sum(log(p)), df=2k, lower.tail=FALSE)` (bhklab/survcomp `R/combine.test.R`; weights ignored). `scipy.stats.combine_pvalues(method='fisher')` (scipy 1.18.1) equals base-R values to < 1e-12 relative on 8 cases, including p = 0 → 0, (1, 0) → 0 and 1e-300 underflow → 0. `false_discovery_control(method='bh')` equals `p.adjust('fdr')`.
- **G3, decided 16/10/2026:** many sites tie at adj p = 0; rule 12 gives the tie-break. R's `order` is stable, so in compat mode ties keep input order.
- **D1–D4, decided 16/10/2026** (decision log DL-018–DL-021): D1 table + name-rule fallback (rule 13); D2 non-missing label, compat `min()` (rule 11); D3 p = 1 + warning (rule 6); D4 `--plots` optional extra.
- Golden p-values still need the shuffled tables (`questions_for_authors.md` Q3; DL-023).
