#!/bin/bash
# Builds the end-counts golden fixture (DL-038): subsample of nanoPARE's test set flower_1 (5P, BODY) + outputs of the unchanged original on it.
# Needs: the nanoPARE clone with the 22/10 test run (temp/flower_1.{5P,BODY}/flower_1.sam), python3, samtools, bedtools.
set -euo pipefail
NP=${NP:-/home/nv358/rds/nanoPARE}
PY=$NP/scripts/python_scripts
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$(mktemp -d)}
cd $HERE
mkdir -p input golden

cp $NP/resources/genome.fasta input/genome.fasta
samtools faidx input/genome.fasta

subsample() { # library fraction
    python3 -I $HERE/subsample_sam.py $NP/temp/flower_1.$1/flower_1.sam $2 > $WORK/flower_1.$1.sam
    samtools view -b --no-PG -o input/flower_1.$1.bam $WORK/flower_1.$1.sam
}
subsample 5P 0.30
subsample BODY 0.15

# unchanged readmapIO.py with the published non-bias flags (endMap.sh l.495-513) + post-processing (l.521-533)
cd $WORK
s=flower_1
python3 $PY/readmapIO.py -S $s -I $WORK/$s.5P.sam -R 5P -F $HERE/input/genome.fasta --softclip_type 5p --untemp_out G --secondary --allow_naive > $HERE/golden/readmapIO_5P.log 2>&1
python3 $PY/readmapIO.py -S $s -I $WORK/$s.BODY.sam -R BODY -F $HERE/input/genome.fasta --secondary --allow_nonstranded --allow_naive > $HERE/golden/readmapIO_BODY.log 2>&1
for lt in 5P BODY; do
    bedtools sort -i ${s}_${lt}_plus.bedgraph > $HERE/golden/$s.$lt.plus.bedgraph
    bedtools sort -i ${s}_${lt}_minus.bedgraph > $HERE/golden/$s.$lt.minus.bedgraph
done
bedtools sort -i ${s}_5P_plus_untemp.bedgraph | awk '{ printf $1"\t"$2"\t"$3"\t"$4"\n" }' > $HERE/golden/$s.uG.plus.bedgraph
bedtools sort -i ${s}_5P_minus_untemp.bedgraph | awk '{ printf $1"\t"$2"\t"$3"\t"$4"\n" }' > $HERE/golden/$s.uG.minus.bedgraph
python3 $PY/bedgraph_combine.py -i ${s}_BODY_plus_coverage.bedgraph ${s}_BODY_minus_coverage.bedgraph -o $HERE/golden/$s.BODY.cov.bedgraph

{
    echo "built $(date -Is) on $(hostname)"
    echo "nanoPARE: $(git -C $NP rev-parse --short HEAD) ($(git -C $NP branch --show-current)); readmapIO.py and bedgraph_combine.py identical to 6264d6f: $(git -C $NP diff --quiet 6264d6f -- scripts/python_scripts/readmapIO.py scripts/python_scripts/bedgraph_combine.py && echo yes || echo NO)"
    echo "python3: $(python3 --version 2>&1); $(samtools --version | head -1); $(bedtools --version)"
} > $HERE/golden/provenance.txt
cd $HERE && md5sum input/* golden/*.bedgraph > MD5SUMS
