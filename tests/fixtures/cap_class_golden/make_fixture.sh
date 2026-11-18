#!/bin/bash
# cap-class golden fixture (18/11/2026). Run from this folder.
set -euo pipefail
R=/home/nv358/rds/nanoPARE/results          # golden test-set run of the original (6264d6f + P1-P2)
K=/home/nv358/rds/rds-apis-mellifera-7BbYnVDdZgw/nanopare_nfcore_design/answers/scratch/16_11_2026_cap-class/root_keep/results/EndClass/flower  # rerun keeping intermediates (identical finals)
for s in 1 2 3; do for k in plus minus; do
  gzip -9nc $R/EndGraph/flower_$s/uG.${k}_mask.bedgraph > input/flower_$s.uG.${k}_mask.bedgraph.gz
done; done
cp $R/EndClass/flower/flower.rep.bed input/flower.rep.bed
for f in plus.bedgraph minus.bedgraph uG.plus.bedgraph uG.minus.bedgraph capped.bed noncapped.bed; do
  gzip -9nc $R/EndClass/flower/flower.$f > golden/flower.$f.gz
done
gzip -9nc $K/flower.peaks.bed > golden/flower.peaks.bed.gz   # removed by the original's cleanup; kept in the rerun
find input golden -type f | sort | xargs md5sum > MD5SUMS
