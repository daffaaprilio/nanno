#!/bin/bash

SORT_INDEX_BAM_DIR="sample/sort_bam"
BED_DIR="sample/bed"

MODKIT="/home/daapr/dist_modkit_v0.4.4_251055f/modkit" # install later in matsu

for file in ${SORT_INDEX_BAM_DIR}/*.bam; do
  filename=$(basename "$file" .bam)
  ${MODKIT} pileup $file ${BED_DIR}/${filename}.bed
done