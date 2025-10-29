#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


def load_bed_files():
    """Load BED files for all samples and genes"""

    # Define column names for BED files
    BED_colnames = [
        "chrom", "start_position", "end_position", "base_code", "score", "strand",
        "start_position2", "end_position2", "color", "Nvalid_cov", "percent_modified",
        "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"
    ]

    # Define gene lists
    lpat_genes = ["LPAT1", "LPAT2", "LPAT3", "LPAT4"]
    # dgat_genes = ["DGATa", "DGATb", "DGATc", "DGATd", "DGATe", "DGATf", "DGATg", "DGAT2", "DGATh"]

    # Define sample lists
    samples = ["NIES_2145", "NIES_2146", "NIES_2145_bta1l"]

    # Define Contigs dictionary
    contigs = {
        "LPAT1": "JBEBFO010000001.1",
        "LPAT2": "JBEBFO010000006.1",
        "LPAT3": "JBEBFO010000010.1",
        "LPAT4": "JBEBFO010000010.1"
    }
    #     "DGATa": "JBEBFO010000001.1",
    #     "DGATb": "JBEBFO010000001.1",
    #     "DGATc": "JBEBFO010000018.1",
    #     "DGATd": "JBEBFO010000002.1",
    #     "DGATe": "JBEBFO010000022.1",
    #     "DGATf": "JBEBFO010000004.1",
    #     "DGATg": "JBEBFO010000004.1",
    #     "DGAT2": "JBEBFO010000005.1",
    #     "DGATh": "JBEBFO010000008.1"
    # }

    # Define sample gene BED file directory
    bed_dir = "sample_gene_bed"

    # Dictionary to store all data
    sample_data = {}

    # Loop through each sample and gene, read .bed files
    for sample in samples:
        lpat_dict = {}
        dgat_dict = {}

        for gene in lpat_genes: # + dgat_genes:
            bed_path = os.path.join(bed_dir, f"{sample}_{gene}.bed")
            if os.path.exists(bed_path):
                bed_df = pd.read_csv(bed_path, sep="\t", header=None, names=BED_colnames)

                # Add to LPAT or DGAT dict
                if gene in lpat_genes:
                    lpat_dict[gene] = {"bed_df": bed_df, "contig": contigs[gene]}
                else:
                    dgat_dict[gene] = {"bed_df": bed_df, "contig": contigs[gene]}
            else:
                print(f"Warning: File not found - {bed_path}")

        # Store dictionaries for each sample
        sample_data[sample] = {
            "LPAT": lpat_dict,
            "DGAT": dgat_dict
        }

    return sample_data, samples, lpat_genes, contigs # , dgat_genes


def create_methylation_plots(sample_data, samples, lpat_genes, contigs):# dgat_genes,
    """Create methylation landscape plots for all samples and genes"""

    # Create output directory for visualizations
    output_dir = "visualization_methylation_landscape"
    os.makedirs(output_dir, exist_ok=True)

    for sample in samples:
        for gene in lpat_genes: # + dgat_genes:
            # Get gene data from the appropriate dictionary
            if gene in lpat_genes:
                gene_data = sample_data[sample]["LPAT"].get(gene)
            else:
                gene_data = sample_data[sample]["DGAT"].get(gene)

            if gene_data is None:
                print(f"Warning: No data found for {sample}_{gene}")
                continue

            sample_gene = gene_data["bed_df"]
            contig = gene_data["contig"]

            # Skip if dataframe is empty
            if sample_gene.empty:
                print(f"Warning: Empty dataframe for {sample}_{gene}")
                continue

            # Adjust scale based on the current gene dataframe
            x_axis_length = sample_gene['end_position'].max() - sample_gene['start_position'].min()

            # Let sns.relplot manage the figure creation
            g = sns.relplot(
                data=sample_gene, x='start_position', y='percent_modified',
                row='base_code', row_order=['m', '21839', 'a'],
                hue='base_code', hue_order=['m', 'a', '21839'],
                size='score',
                kind='scatter', height=3, aspect=7,
                legend=True
            )

            # Set title and labels
            g.figure.suptitle(f"{sample}: Methylation in {gene} gene (Contig: {contig})", x=0.5, y=1.02)
            g.set_axis_labels("Base number", "% methylated base")
            g.set(xlim=(sample_gene['start_position'].min(), sample_gene['start_position'].min() + x_axis_length))
            g.set(ylim=(0, 105))
            g._legend.set(title='Depth')

            # Save plot as image
            output_path = os.path.join(output_dir, f"{sample}_{gene}.png")
            plt.savefig(output_path, bbox_inches='tight')
            plt.close()  # Close the figure to free memory

            print(f"Saved plot: {output_path}")


def main():
    """Main function to run the methylation landscape visualization"""

    print("Loading BED files...")
    # sample_data, samples, lpat_genes, dgat_genes, contigs = load_bed_files()
    sample_data, samples, lpat_genes, contigs = load_bed_files()

    print("Creating methylation plots...")
    create_methylation_plots(sample_data, samples, lpat_genes, contigs) # , dgat_genes

    print("Methylation landscape visualization completed!")


if __name__ == "__main__":
    main()
