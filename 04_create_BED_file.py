#!/usr/bin/env python3

import os
import pandas as pd
import glob

# Function to parse BLAST results
def parse_blast_result(filename):
    with open(filename, 'r') as file:
        first_line = file.readline().strip()
    columns = first_line.split('\t')
    chrom = columns[1]
    col9, col10 = int(columns[8]), int(columns[9])
    return chrom, min(col9, col10), max(col9, col10)

# Preparation 1: Load whole genome BED files
BED_colnames = ["chrom", "start_position", "end_position", "base_code", "score", "strand",
                "start_position2", "end_position2", "color", "Nvalid_cov", "percent_modified",
                "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"]

sample_dfs = {
    "NIES_2145": pd.read_csv("sample/bed/Nanoce_2145_sorted.bed", sep='\t', header=None, names=BED_colnames),
    "NIES_2146": pd.read_csv("sample/bed/Nanoce_2146_sorted.bed", sep='\t', header=None, names=BED_colnames),
    "NIES_2145_bta1l": pd.read_csv("sample/bed/Nanoce_2145_bta1l_sorted.bed", sep='\t', header=None, names=BED_colnames)
}

# Preparation 2: Automatically load coordinates from ALL BLAST result files
gene_coords = {}
blast_files = glob.glob("*_blast.txt")  # Find all files ending with _blast.txt

for blast_file in blast_files:
    gene_name = blast_file.replace("_blast.txt", "")  # Extract gene name from filename
    chrom, start, end = parse_blast_result(blast_file)
    gene_coords[gene_name] = {"chrom": chrom, "start": start, "end": end}
    print(f"Loaded {gene_name}: {chrom}:{start}-{end}")

# Create output folder
output_dir = "sample_gene_bed"
os.makedirs(output_dir, exist_ok=True)

# Filter and save for each sample and gene
for sample_name, df in sample_dfs.items():
    for gene_name, coords in gene_coords.items():
        filtered_df = df[
            (df['chrom'] == coords['chrom']) &
            (df['start_position'] >= coords['start']) &
            (df['end_position'] <= coords['end'])
        ].reset_index(drop=True)

        # Save as variable and export as BED
        var_name = f"{sample_name}_{gene_name}"
        globals()[var_name] = filtered_df

        bed_path = f"{output_dir}/{var_name}.bed"
        filtered_df.to_csv(bed_path, sep="\t", index=False, header=False)

print(f"Processed {len(gene_coords)} genes for {len(sample_dfs)} samples")
