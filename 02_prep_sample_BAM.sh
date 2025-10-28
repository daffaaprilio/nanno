#!/bin/bash

BAM_DIR="sample/bam"
SORT_INDEX_BAM_DIR="sample/sort_bam"

SAMTOOLS="/home/daffa/bin/samtools" # "/home/daapr/samtools/bin/samtools"
IGVTOOLS="/home/daffa/IGV_2.17.4/igvtools" # install later in matsu (when in lab)

for file in "${BAM_DIR}"/*.bam; do
  filename=$(basename "$file" .bam)
  ${SAMTOOLS} sort -o ${SORT_INDEX_BAM_DIR}/${filename}_sorted.bam $file && ${SAMTOOLS} index ${SORT_INDEX_BAM_DIR}/${filename}_sorted.bam
done

for file in ${SORT_INDEX_BAM_DIR}/*.bam; do
  filename=$(basename "$file" .bam)
  
  ${IGVTOOLS} index ${file}
done

