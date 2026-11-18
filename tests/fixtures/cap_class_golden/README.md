# cap-class golden fixture

Golden outputs of nanoPARE's EndClass (numeric part) on nanoPARE's test set, group `flower` (3 replicates), for `cap-class --nanopare_compat` and the default-vs-golden check.

## Contents
| Path | What |
|---|---|
| `input/flower_{1,2,3}.uG.{plus,minus}_mask.bedgraph.gz` | EndGraph cap tracks (original format: G one nt upstream of the 5′ end) |
| `input/flower.rep.bed` | the original's replicable regions (`flower.rep.bed`) |
| `golden/flower.{plus,minus}.bedgraph.gz`, `golden/flower.uG.{plus,minus}.bedgraph.gz` | the original's replicate sums |
| `golden/flower.peaks.bed.gz` | trimmed regions + peak (deleted by the original's cleanup; taken from the rerun that kept intermediates, whose final outputs are identical to golden) |
| `golden/flower.{capped,noncapped}.bed.gz` | the original's split, before renaming |
| `make_fixture.sh`, `MD5SUMS` | how the files were made (18/11/2026) and their checksums |

5P inputs are shared with `../scale_factor_test_set/input/` (`flower_*.5P.{plus,minus}_mask.bedgraph.gz`, md5-identical to the original's EndGraph files); lengths from `../kde_peaks_golden/length.table`.

## Source
`results/EndClass/flower/` and `results/EndGraph/flower_*/` of the original's test-set run (clone `~/rds/nanoPARE`, 6264d6f + P1–P2); EndClass rerun 16/11 reproduced 14/14 golden files.

## Expected
- Compat: the 4 sum bedgraphs byte-identical; `peaks.bed` columns 1–6 identical (sorted); `capped.bed` / `noncapped.bed` identical in columns 1–3, 6, 10 (sorted); 323 / 279.
- Default: inputs are the original's uG tracks moved to the 5′ end (+ strand +1, − strand −1) and capped at the 5P value per nt (approximates one count per read); see `tests/test_cap_class_golden.py`.
