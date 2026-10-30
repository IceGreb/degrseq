# cleavage-stats golden fixture

Golden test data for `cleavage-stats --nanopare_compat`.

## Contents
| Path | What |
|---|---|
| `input/<sample>.real.w{20,50}.tsv` | EndCut step-1 tables of the real sRNA set (flower_1–3) |
| `input/<sample>.shuffled.w{20,50}.tsv.gz` | the 1,000 shuffled step-1 tables pooled in set order, one header (DL-050) |
| `input/srna_table.tsv` | `anno.miRNA.tsv` + `anno.tasiRNA.tsv` names with their class |
| `input/precursors.tsv`, `input/validated_targets.tsv` | `miRBase21_and_TAIR10_miRNA_tasiRNA_precursors.txt`, `sRNA_targets_5pRACE.tsv`, unchanged |
| `input/gene_aliases.tsv` | `gene_aliases_20130130.txt`, AT2G rows only (file order kept, so the first row per gene is unchanged) |
| `golden/<sample>.{miRNA,tasiRNA}.detected_sites.tsv` | `<sample>.anno.{mir,tas}.allen.detected.sites.tsv` of the original |
| `MD5SUMS` | checksums |

## How it was made (`make_fixture.sh`, 30/10/2026)
From the golden run of the original on its test set (`nanopare_golden_run/EndCut/`, 29/10/2026):
EndCut step 1 with `drivers/03_step1.sh`, then `EndCut_step2.v03.git.R` from nanoPARE `8ce221c`
(6264d6f + P1–P3; P3 = `method="radix"` at l.79) with `drivers/04_step2.sh` (R 4.5.3, survcomp 1.60.1, `LC_ALL=C`).

`adj.p.val` can differ from R in the 15th significant digit (scipy's chi-square tail vs R's
`pchisq`; 5 of 41 values). The test compares it with a relative tolerance of 1e-12 and every
other field as text (DL-052).
