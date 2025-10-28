#!/bin/bash

RAW_DIR="raw/bam.d"
OUT_DIR="sample/bam"

SAMTOOLS="/home/daffa/bin/samtools" # "/home/daapr/samtools/bin/samtools"
# igvtools="/home/daapr/IGV_2.17.4/igvtools"
# modkit="/home/daapr/dist_modkit_v0.4.4_251055f/modkit"

${SAMTOOLS} merge -o ${OUT_DIR}/Nanoce_2145.bam ${RAW_DIR}/r0039n.bam ${RAW_DIR}/r0041n.bam ${RAW_DIR}/r0043n.bam
${SAMTOOLS} merge -o ${OUT_DIR}/Nanoce_2145_bta1l.bam ${RAW_DIR}/r0037n.bam ${RAW_DIR}/r0038n.bam
${SAMTOOLS} merge -o ${OUT_DIR}/Nanoce_2146.bam ${RAW_DIR}/r0040n.bam ${RAW_DIR}/r0042n.bam

