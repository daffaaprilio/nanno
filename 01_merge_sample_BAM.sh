#!/bin/bash

raw_bam_dir="raw/bam.d"
sample_bam_dir="nanno/bam"

SAMTOOLS="/home/daffa/bin/samtools" # "/home/daapr/samtools/bin/samtools"
# igvtools="/home/daapr/IGV_2.17.4/igvtools"
# modkit="/home/daapr/dist_modkit_v0.4.4_251055f/modkit"

${SAMTOOLS} merge -o nanno/bam/Nanoce_2145.bam raw/bam.d/r0039n.bam raw/bam.d/r0041n.bam raw/bam.d/r0043n.bam
${SAMTOOLS} merge -o nanno/bam/Nanoce_2145_bta1l.bam raw/bam.d/r0037n.bam raw/bam.d/r0038n.bam
${SAMTOOLS} merge -o nanno/bam/Nanoce_2146.bam raw/bam.d/r0040n.bam raw/bam.d/r0042n.bam

