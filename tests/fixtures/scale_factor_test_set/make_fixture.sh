#!/bin/bash
# Builds the scale-factor test-set fixture: EndGraph inputs of the 22/10 nanoPARE test run (3 samples) + gene regions.
# Needs: the nanoPARE clone with results/EndGraph/flower_{1,2,3}, python3, gzip.
set -euo pipefail
NP=${NP:-/home/nv358/rds/nanoPARE}
HERE=$(cd "$(dirname "$0")" && pwd)
cd $HERE
mkdir -p input
for s in flower_1 flower_2 flower_3; do
    for f in 5P.plus_mask 5P.minus_mask BODY.plus BODY.minus; do
        gzip -n -9 -c $NP/results/EndGraph/$s/$f.bedgraph > input/$s.$f.bedgraph.gz
    done
done
python3 -I make_regions.py $NP/resources/annotation.gff | gzip -n -9 > input/regions.bed.gz
md5sum input/* > MD5SUMS
