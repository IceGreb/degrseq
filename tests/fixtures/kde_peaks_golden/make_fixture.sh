#!/bin/bash
# Golden outputs of nanoPARE's EndGraph on its test set (clone ~/rds/nanoPARE, branch golden-p1, 22/10 run); inputs are shared with ../scale_factor_test_set/input
set -euo pipefail
G=${1:?usage: make_fixture.sh <nanoPARE>/results/EndGraph}
cd "$(dirname "$0")"
for s in flower_1 flower_2 flower_3; do
  gzip -n -c $G/$s/$s.end_features.bed > golden/$s.end_features.bed.gz
  for st in plus minus; do gzip -n -c $G/$s/5P.${st}_smooth.bedgraph > golden/$s.5P.${st}_smooth.bedgraph.gz; done
done
md5sum golden/* length.table scale_factors.tsv > MD5SUMS
