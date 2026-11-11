# scale-factor test-set fixture

Inputs for checking `scale-factor` (default mode; no compat, DL-058) against the audit's measured values on nanoPARE's test set.

## Contents
| Path | What |
|---|---|
| `input/flower_{1,2,3}.5P.{plus,minus}_mask.bedgraph.gz` | 5P ends after the TSO mask, as EndGraph passes them to `gtf_quantify.py` |
| `input/flower_{1,2,3}.BODY.{plus,minus}.bedgraph.gz` | BODY ends (runs of equal value merged, as the original writes them) |
| `input/regions.bed.gz` | gene regions: exons of each transcript ∪ 100 nt before/after, merged per gene and strand; gene = transcript ID before `.` |
| `make_regions.py`, `make_fixture.sh` | how the inputs were made (11/11/2026) |
| `MD5SUMS` | checksums of `input/` |

## Source
`results/EndGraph/flower_{1,2,3}/` of the 22/10 test run of the original (clone `~/rds/nanoPARE`, branch `golden-p1`, commit `35f4bed` = 6264d6f + two patches) and `resources/annotation.gff`. `make_regions.py` stands in for the pipeline's gffread + awk + `bedtools merge` step.

## Expected (`tests/test_scale_factor_test_set.py`)
| Sample | factor | E | B | golden (original) |
|---|---|---|---|---|
| flower_1 | 0.3871 | 138,729 | 160,859 | 0.3924 |
| flower_2 | 0.1903 | 214,506 | 121,465 | 0.2100 |
| flower_3 | 0.2103 | 219,022 | 139,615 | 0.2121 |

Differences from golden: known issues 43–45 and gene length = region size (spec quirk table). flower_1's B includes AT2G40740 (2 reads), which the original drops (#43).
