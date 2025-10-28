#!/bin/bash

bam_dir="nanno/bam"
sort_index_bam_dir="nanno/sort_bam"

samtools="/home/daapr/samtools/bin/samtools"
igvtools="/home/daapr/IGV_2.17.4/igvtools"

for file in "$bam_dir"/*.bam; do
  filename=$(basename "$file" .bam)
  "$samtools" sort -o "$sort_index_bam_dir"/${filename}_sorted.bam $file && "$samtools" index "$sort_index_bam_dir"/${filename}_sorted.bam
done

for file in "$sort_index_bam_dir"/*.bam; do
  filename=$(basename "$file" .bam)
  
  "$igvtools" index "$file"
done

