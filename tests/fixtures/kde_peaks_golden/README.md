# kde-peaks golden fixture

Golden outputs of nanoPARE's EndGraph steps 2b–3 on nanoPARE's test set, for `kde-peaks --nanopare_compat` (byte comparison) and the default-vs-golden check.

## Contents
| Path | What |
|---|---|
| `golden/flower_{1,2,3}.end_features.bed.gz` | the original's `<sample>.end_features.bed` |
| `golden/flower_{1,2,3}.5P.{plus,minus}_smooth.bedgraph.gz` | the original's smoothed density (`--write_density` output) |
| `length.table` | the test set's chromosome lengths (`resources/length.table`) |
| `scale_factors.tsv` | the golden scale factors the original used (0.3924 / 0.21 / 0.2121) |
| `make_fixture.sh`, `MD5SUMS` | how the files were made (16/11/2026) and their checksums |

Inputs (5P masked ends, BODY ends) are shared with `../scale_factor_test_set/input/` and are identical to the original's EndGraph inputs (md5 checked 16/11/2026).

## Source
`results/EndGraph/flower_{1,2,3}/` of the 22/10/2026 test-set run of the original (clone `~/rds/nanoPARE`, branch `golden-p1`, 6264d6f + P1–P2; the EndGraph helpers are unchanged by the patches). The audit (`nanopare_nfcore_design/answers/12_11_2026_kde-peaks_audit.md` §7) re-ran the original on these inputs and reproduced these files byte for byte.

## Expected
- Compat: density files byte-identical; `end_features.bed` byte-identical after sorting by chromosome, start and strand (spec G2: flower_2 has one pair of rows with equal start whose order the original leaves undefined).
- Default: see `tests/test_kde_peaks_golden.py` (counts from the spec's quirk table).
