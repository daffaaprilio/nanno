#!/bin/bash

sort_index_bam_dir="nanno/sort_bam"
bed_dir="nanno/bed"

modkit="/home/daapr/dist_modkit_v0.4.4_251055f/modkit"

for file in "$sort_index_bam_dir"/*.bam; do
  filename=$(basename "$file" .bam)
  "$modkit" pileup $file "$bed_dir"/${filename}.bed
done