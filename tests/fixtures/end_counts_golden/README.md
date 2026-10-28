# end-counts golden fixture

Golden test data for `end-counts --nanopare_compat` (decision DL-038).

## Contents
| Path | What |
|---|---|
| `input/genome.fasta` (+ `.fai`) | nanoPARE test genome: 1 Mb of Arabidopsis chr2 (`Ath_chr2`) |
| `input/flower_1.5P.bam` | 5P subsample: 53,928 reads (57,750 records) |
| `input/flower_1.BODY.bam` | BODY subsample: 54,990 reads (62,887 records) |
| `golden/flower_1.{5P,BODY}.{plus,minus}.bedgraph` | 5′ ends, as `endMap.sh` writes them |
| `golden/flower_1.uG.{plus,minus}.bedgraph` | 5P cap track (untemplated G) |
| `golden/flower_1.BODY.cov.bedgraph` | BODY coverage, both strands summed |
| `golden/readmapIO_{5P,BODY}.log`, `golden/provenance.txt` | run logs and tool versions |
| `MD5SUMS` | checksums of `input/` and `golden/` |

## How it was made (`make_fixture.sh`, 28/10/2026)
1. **Source reads:** STAR output of nanoPARE's test set, sample flower_1. These are the exact SAMs `readmapIO.py` read in the 22/10 test run of the original (`temp/flower_1.{5P,BODY}/flower_1.sam` in the nanoPARE clone).
2. **Subsample (`subsample_sam.py`):** keeps every multimapped read (NH > 1) with all its alignments, plus unique reads whose name hash (MD5) falls below 0.30 (5P) or 0.15 (BODY). Deterministic: the same reads are kept on every run.
3. **Golden output:** the unchanged original (`readmapIO.py`, `bedgraph_combine.py` identical to nanoPARE 6264d6f), with the published flags (`endMap.sh` l.495–513) and the same post-processing (l.521–533), run **on the subsample**. Multimapper sharing depends on the whole input's coverage, so the full-sample outputs can't serve as the reference.

Rebuilding gives byte-identical output (checked with a different `PYTHONHASHSEED`). Python ≥ 3.6 is needed, because the original's tie order follows dict order.
