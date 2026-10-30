#!/bin/bash
# Builds this fixture from the golden run of the original (nanopare_golden_run/EndCut, 29/10/2026).
set -euo pipefail
G=${1:-/home/nv358/rds/rds-apis-mellifera-7BbYnVDdZgw/nanopare_golden_run/EndCut}
cd "$(dirname "$0")"
export LC_ALL=C
R=$G/results/transcript_bedgraph_capmasked
A=$G/Ath_annotations
rm -rf input golden && mkdir -p input golden

# sRNA table: the two class tables of the original, names whitespace-stripped
{ printf 'sRNA\tclass\n'
  tail -n +2 $A/anno.miRNA.tsv | awk -F'\t' -v OFS='\t' '{gsub(/ +$/, "", $1); print $1, "miRNA"}'
  tail -n +2 $A/anno.tasiRNA.tsv | awk -F'\t' -v OFS='\t' '{gsub(/ +$/, "", $1); print $1, "tasiRNA"}'
} > input/srna_table.tsv
cp $A/miRBase21_and_TAIR10_miRNA_tasiRNA_precursors.txt input/precursors.tsv
cp $A/sRNA_targets_5pRACE.tsv input/validated_targets.tsv
# aliases: chr2 genes only (the test transcriptome), file order kept so the first row per gene is unchanged
awk -F'\t' 'NR == 1 || $1 ~ /^AT2G/' $A/gene_aliases_20130130.txt > input/gene_aliases.tsv

for s in flower_1 flower_2 flower_3; do
    for w in 20 50; do
        cp $R/$s/anno.mir.tas.fa.GSTAr_${s}.transcript.capmasked_win_${w}_detectedSites.txt input/$s.real.w$w.tsv
        # all 1,000 shuffled tables pooled in set order, one header (DL-050)
        { head -1 input/$s.real.w$w.tsv
          for j in $(seq 0 999); do
              tail -n +2 $R/$s/anno.mir.tas.$j.shuffled.fa.GSTAr_${s}.transcript.capmasked_win_${w}_detectedSites.txt
          done
        } | gzip -n > input/$s.shuffled.w$w.tsv.gz
    done
    cp $R/$s/$s.anno.mir.allen.detected.sites.tsv golden/$s.miRNA.detected_sites.tsv
    cp $R/$s/$s.anno.tas.allen.detected.sites.tsv golden/$s.tasiRNA.detected_sites.tsv
done
md5sum input/* golden/* > MD5SUMS
