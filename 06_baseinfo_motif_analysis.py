#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import subprocess
import sys


def setup_data():
    """Setup gene lists, samples, and contigs configuration"""

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
        "LPAT4": "JBEBFO010000010.1",
        "DGATa": "JBEBFO010000001.1",
        "DGATb": "JBEBFO010000001.1",
        "DGATc": "JBEBFO010000018.1",
        "DGATd": "JBEBFO010000002.1",
        "DGATe": "JBEBFO010000022.1",
        "DGATf": "JBEBFO010000004.1",
        "DGATg": "JBEBFO010000004.1",
        "DGAT2": "JBEBFO010000005.1",
        "DGATh": "JBEBFO010000008.1"
    }

    return lpat_genes, samples, contigs
    # return lpat_genes, dgat_genes, samples, contigs


def load_bed_files(samples, lpat_genes, contigs):
# def load_bed_files(samples, lpat_genes, dgat_genes, contigs):
    """Load BED files for all samples and genes"""

    BED_colnames = [
        "chrom", "start_position", "end_position", "base_code", "score", "strand",
        "start_position2", "end_position2", "color", "Nvalid_cov", "percent_modified",
        "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"
    ]

    bed_dir = "sample_gene_bed"
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

        sample_data[sample] = {
            "LPAT": lpat_dict,
            # "DGAT": dgat_dict
        }

    return sample_data, BED_colnames


def generate_sense_base_info(samples, lpat_genes):
# def generate_sense_base_info(samples, lpat_genes, dgat_genes):
    """Generate sense base information using bedtools getfasta"""

    reference_genome = "ref/Nanoce_C018.fna"
    bed_dir = "sample_gene_bed"

    print("Generating sense base information...")

    for sample in samples:
        for gene in lpat_genes: # + dgat_genes:
            bed_path = f"{bed_dir}/{sample}_{gene}.bed"
            sensebase_bed_path = f"{bed_dir}/sensebase_{sample}_{gene}.bed"

            if os.path.exists(bed_path):
                try:
                    # Run bedtools getfasta
                    cmd = f"bedtools getfasta -fi {reference_genome} -bed {bed_path} -bedOut"
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                    with open(sensebase_bed_path, 'w') as f:
                        f.write(result.stdout)

                    print(f"Generated: {sensebase_bed_path}")
                except Exception as e:
                    print(f"Error generating {sensebase_bed_path}: {e}")


def antisense(baseinfo_df):
    """Convert sense bases to antisense bases"""
    # Complementary base mapping
    complement_map = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                      'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    # Apply map directly on the column (vectorized, much faster)
    baseinfo_df["antisense_base"] = baseinfo_df["sense_base"].map(complement_map)

    return baseinfo_df


def combine_base_information(samples, lpat_genes): #, dgat_genes):
    """Combine sense and antisense base information"""

    BED_sensebase_colnames = [
        "chrom", "start_position", "end_position", "base_code", "score", "strand",
        "start_position2", "end_position2", "color", "Nvalid_cov", "percent_modified",
        "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall", "sense_base"
    ]

    bed_dir = "sample_gene_bed"
    baseinfo_data = {}

    print("Combining base information...")

    for sample in samples:
        baseinfo_data[sample] = {}
        for gene in lpat_genes: # + dgat_genes:
            sensebase_bed_path = f"{bed_dir}/sensebase_{sample}_{gene}.bed"

            if os.path.exists(sensebase_bed_path):
                df = pd.read_csv(sensebase_bed_path, sep='\t', header=None, names=BED_sensebase_colnames)
                df_antisense = antisense(df)
                baseinfo_data[sample][gene] = df_antisense
                print(f"Processed: {sample}_{gene}_baseinfo")

    return baseinfo_data


def find_methylation_motif(df_baseinfo):
    """Identify methylation motifs (CpG, CHG, CHH)"""

    df_motifinfo = df_baseinfo.copy()
    # Remove duplicate C, which is 4mC
    df_motifinfo = df_motifinfo[df_motifinfo['base_code'] != '21839']
    # Remove duplicate base which has 6mA
    duplicate_base_1 = df_motifinfo['start_position'] == df_motifinfo['start_position'].shift(1)
    duplicate_base_min1 = df_motifinfo['start_position'] == df_motifinfo['start_position'].shift(-1)
    remove_duplicate_base = duplicate_base_1 & (df_motifinfo['base_code'] == 'a')
    remove_duplicate_base = duplicate_base_min1 & (df_motifinfo['base_code'] == 'a')
    df_motifinfo = df_motifinfo[~remove_duplicate_base].reset_index(drop=True)

    # Create 'sense_motif' and 'antisense_motif' columns
    def classify_motif(base_seq, next_base1, next_base2):
        """Classifies the methylation motif based on base sequences."""
        if base_seq in ['C','c']:
            if next_base1 in ['G','g']:
                return 'CpG'
            elif next_base1 in ['A', 'a', 'T', 't', 'C', 'c']:
                if next_base2 in ['G', 'g']:
                    return 'CHG'
                elif next_base2 in ['A', 'a', 'T', 't', 'C', 'c']:
                    return 'CHH'
        return np.nan

    # Shift columns to access the next base(s)
    df_motifinfo['next_sense_base1'] = df_motifinfo['sense_base'].shift(-1)
    df_motifinfo['next_sense_base2'] = df_motifinfo['sense_base'].shift(-2)

    df_motifinfo['sense_motif'] = df_motifinfo.apply(
        lambda row: classify_motif(row['sense_base'], row['next_sense_base1'], row['next_sense_base2']), axis=1
    )

    df_motifinfo['next_antisense_base1'] = df_motifinfo['antisense_base'].shift(-1)
    df_motifinfo['next_antisense_base2'] = df_motifinfo['antisense_base'].shift(-2)

    df_motifinfo['antisense_motif'] = df_motifinfo.apply(
        lambda row: classify_motif(row['antisense_base'], row['next_antisense_base1'], row['next_antisense_base2']), axis=1
    )

    # Drop helper columns used for shifting
    df_motifinfo.drop(['next_sense_base1', 'next_sense_base2', 'next_antisense_base1', 'next_antisense_base2'], axis=1, inplace=True)

    return df_motifinfo


def process_motif_analysis(baseinfo_data, samples, lpat_genes): #, dgat_genes):
    """Execute motif analysis and filter by motif types"""

    motif_data = {}
    motifs = ["CpG", "CHG", "CHH"]

    print("Processing motif analysis...")

    for sample in samples:
        motif_data[sample] = {}
        for gene in lpat_genes: # + dgat_genes:
            if gene in baseinfo_data[sample]:
                # Execute find motif info function
                motifinfo = find_methylation_motif(baseinfo_data[sample][gene])
                motif_data[sample][gene] = {"motifinfo": motifinfo}

                # Filter by motif types
                for motif in motifs:
                    subset_df = motifinfo.loc[
                        ((motifinfo["sense_motif"] == motif) | (motifinfo["antisense_motif"] == motif)) &
                        (motifinfo["base_code"] == "m")
                    ]
                    motif_data[sample][gene][motif] = subset_df

    return motif_data


def create_methylation_landscape_plots(motif_data, samples, lpat_genes, contigs): # dgat_genes, contigs):
    """Create methylation landscape plots per motif"""

    motifs = ["CpG", "CHG", "CHH"]
    output_dir = "visualization_motif_analysis/methylation_landscape"
    os.makedirs(output_dir, exist_ok=True)

    print("Creating methylation landscape plots...")

    for sample in samples:
        for gene in lpat_genes: # + dgat_genes:
            if gene in motif_data[sample]:
                for motif in motifs:
                    if motif in motif_data[sample][gene]:
                        contig = contigs[gene]
                        sample_gene_motif = motif_data[sample][gene][motif]

                        if sample_gene_motif.empty:
                            continue

                        # Adjust scale based on the current gene dataframe
                        x_axis_length = sample_gene_motif['end_position'].max() - sample_gene_motif['start_position'].min()

                        # Let sns.relplot manage the figure creation
                        g = sns.relplot(
                            data=sample_gene_motif, x='start_position', y='percent_modified', size='score',
                            kind='scatter', height=3, aspect=3.5,
                            legend=True
                        )

                        # Set title and labels
                        g.figure.suptitle(f"{sample}: {motif} Methylation in {gene} gene (Contig: {contig})", x=0.5, y=1.02)
                        g.set_axis_labels("Base number", "% methylated base")
                        g.set(xlim=(sample_gene_motif['start_position'].min(), sample_gene_motif['start_position'].min() + x_axis_length))
                        g.set(ylim=(0,105))
                        g._legend.set(title='Depth')

                        # Save plot as image
                        output_path = f"{output_dir}/{sample}_{gene}_{motif}.png"
                        plt.savefig(output_path, bbox_inches='tight')
                        plt.close()
                        print(f"Saved: {output_path}")


def generate_whole_genome_base_info():
    """Generate sense base information for whole genome BED files"""

    whole_genome_bed_files = [
        "Nanoce_2145_sorted.bed",
        "Nanoce_2146_sorted.bed",
        "Nanoce_2145_bta1l_sorted.bed"
    ]

    reference_genome = "ref/Nanoce_C018.fna"

    print("Generating whole genome base information...")

    for bed_file in whole_genome_bed_files:
        input_bed = f"sample/bed/{bed_file}"
        output_bed = f"sample/bed/sensebase_{bed_file}"

        if os.path.exists(input_bed):
            try:
                cmd = f"bedtools getfasta -fi {reference_genome} -bed {input_bed} -bedOut"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                with open(output_bed, 'w') as f:
                    f.write(result.stdout)

                print(f"Generated: {output_bed}")
            except Exception as e:
                print(f"Error generating {output_bed}: {e}")


def load_whole_genome_base_info():
    """Load whole genome BED files with base information"""

    BED_sensebase_colnames = [
        "chrom", "start_position", "end_position", "base_code", "score", "strand",
        "start_position2", "end_position2", "color", "Nvalid_cov", "percent_modified",
        "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall", "sense_base"
    ]

    samples_wholegenome = [
        'Nanoce_2145',
        'Nanoce_2146',
        'Nanoce_2145_bta1l'
    ]

    wholegenome_data = {}

    print("Loading whole genome base information...")

    for sample in samples_wholegenome:
        sensebase_bed_path = f"sample/bed/sensebase_{sample}_sorted.bed"

        if os.path.exists(sensebase_bed_path):
            df = pd.read_csv(sensebase_bed_path, sep='\t', header=None, names=BED_sensebase_colnames)
            df_antisense = antisense(df)
            wholegenome_data[sample] = df_antisense
            print(f"Loaded: {sample}_baseinfo")

    return wholegenome_data


def find_CN_dimer(df_baseinfo):
    """Identify CN type (for 5mC only)"""

    df_CNinfo = df_baseinfo.copy()
    # Remove duplicate C, which is 4mC
    df_CNinfo = df_CNinfo[df_CNinfo['base_code'] != '21839']
    # Remove duplicate base which has 6mA
    duplicate_base_1 = df_CNinfo['start_position'] == df_CNinfo['start_position'].shift(1)
    duplicate_base_min1 = df_CNinfo['start_position'] == df_CNinfo['start_position'].shift(-1)
    is_6mA = df_CNinfo['base_code'].astype(str) == 'a'
    # apply the filter
    remove_duplicate_base = (duplicate_base_1 | duplicate_base_min1) & is_6mA

    df_CNinfo = df_CNinfo[~remove_duplicate_base].reset_index(drop=True)

    # Ensure sorting by chrom and start_position
    df_CNinfo.sort_values(['chrom', 'start_position'], inplace=True)

    # Get next bases efficiently using vectorized `groupby().shift()`
    df_CNinfo['next_sense_base1'] = df_CNinfo.groupby('chrom')['sense_base'].shift(-1).str.upper()
    df_CNinfo['next_antisense_base1'] = df_CNinfo.groupby('chrom')['antisense_base'].shift(1).str.upper()

    # Define base pair mapping
    base_map = {
        'G': 'CG', 'C': 'CC',
        'A': 'CA', 'T': 'CT'
    }

    # Vectorized motif classification
    df_CNinfo['sense_motif'] = np.where(
        df_CNinfo['sense_base'].str.upper() == 'C',
        df_CNinfo['next_sense_base1'].map(base_map),
        np.nan
    )

    df_CNinfo['antisense_motif'] = np.where(
        df_CNinfo['antisense_base'].str.upper() == 'C',
        df_CNinfo['next_antisense_base1'].map(base_map),
        np.nan
    )

    # Drop temporary columns
    df_CNinfo.drop(columns=['next_sense_base1', 'next_antisense_base1'], inplace=True)

    return df_CNinfo


def calculate_CN(df):
    """Calculate CN dimer ratios"""

    df_filtered = df[df['sense_motif'].notna() | df['antisense_motif'].notna()]

    # Create a long-form dataframe for both strands
    dimer_counts = df_filtered.melt(
        value_vars=['sense_motif', 'antisense_motif'],
        var_name='Strand',
        value_name='Dimer_type'
    ).dropna()

    # Map the original rows to the melted dataframe to bring back Nvalid_cov and Nmod
    dimer_counts['Nvalid_cov'] = dimer_counts.index.map(df['Nvalid_cov'])
    dimer_counts['Nmod'] = dimer_counts.index.map(df['Nmod'])

    # Group by Dimer_type and calculate Frequency and methylated_Frequency
    result = dimer_counts.groupby('Dimer_type').agg(
        Frequency=('Nvalid_cov', 'sum'),
        methylated_Frequency=('Nmod', 'sum')
    ).reset_index()

    # calculate ratio
    result['ratio'] = result['methylated_Frequency']/(result['methylated_Frequency']+result['Frequency'])

    return result


def process_CN_analysis(wholegenome_data):
    """Process CN dimer analysis and create visualizations"""

    samples_wholegenome = list(wholegenome_data.keys())
    dimer_data = {}

    print("Processing CN dimer analysis...")

    for sample in samples_wholegenome:
        CN_info = find_CN_dimer(wholegenome_data[sample])
        dimer_info = calculate_CN(CN_info)
        dimer_data[sample] = dimer_info

    # Create visualizations
    output_dir = "visualization_motif_analysis/CN_analysis"
    os.makedirs(output_dir, exist_ok=True)

    for sample in samples_wholegenome:
        plt.figure(figsize=(7,4))

        sns.scatterplot(data=dimer_data[sample],
                        x='Frequency',
                        y='ratio',
                        zorder=2)

        for i in range(dimer_data[sample].shape[0]):
            plt.text(x=dimer_data[sample]['Frequency'].iloc[i],
                    y=dimer_data[sample]['ratio'].iloc[i],
                    s=dimer_data[sample]['Dimer_type'].iloc[i])

        plt.ylabel("5mCN to total CN ratio")
        plt.xlabel("Total occurrence of CN")
        plt.title(f"{sample}")
        plt.grid(axis='both', zorder=0, color='lightgrey')

        sns.despine()

        output_path = f"{output_dir}/CN_dimer_{sample}.png"
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f"Saved CN analysis: {output_path}")

    return dimer_data


def find_CNN_trimer(df_baseinfo):
    """Identify CNN type (the next 2 bases after methylated)"""

    df_CNNinfo = df_baseinfo.copy()
    # Remove duplicate C, which is 4mC
    df_CNNinfo = df_CNNinfo[df_CNNinfo['base_code'] != '21839']
    # Remove duplicate base which has 6mA
    duplicate_base_1 = df_CNNinfo['start_position'] == df_CNNinfo['start_position'].shift(1)
    duplicate_base_min1 = df_CNNinfo['start_position'] == df_CNNinfo['start_position'].shift(-1)
    is_6mA = df_CNNinfo['base_code'].astype(str) == 'a'
    # apply the filter
    remove_duplicate_base = (duplicate_base_1 | duplicate_base_min1) & is_6mA

    df_CNNinfo = df_CNNinfo[~remove_duplicate_base].reset_index(drop=True)

    # Ensure sorting by chrom and start_position
    df_CNNinfo.sort_values(['chrom', 'start_position'], inplace=True)

    # Get next two bases efficiently using vectorized `groupby().shift()`
    df_CNNinfo['next_sense_base1'] = df_CNNinfo.groupby('chrom')['sense_base'].shift(-1).str.upper()
    df_CNNinfo['next_sense_base2'] = df_CNNinfo.groupby('chrom')['sense_base'].shift(-2).str.upper()

    df_CNNinfo['next_antisense_base1'] = df_CNNinfo.groupby('chrom')['antisense_base'].shift(1).str.upper()
    df_CNNinfo['next_antisense_base2'] = df_CNNinfo.groupby('chrom')['antisense_base'].shift(2).str.upper()

    # Define CNN mapping
    base_map = {
        'GG': 'CGG', 'GC': 'CGC', 'GA': 'CGA', 'GT': 'CGT',
        'CG': 'CCG', 'CC': 'CCC', 'CA': 'CCA', 'CT': 'CCT',
        'AG': 'CAG', 'AC': 'CAC', 'AA': 'CAA', 'AT': 'CAT',
        'TG': 'CTG', 'TC': 'CTC', 'TA': 'CTA', 'TT': 'CTT',
    }

    # Vectorized CNN motif classification
    df_CNNinfo['sense_motif'] = np.where(
        df_CNNinfo['sense_base'].str.upper() == 'C',
        (df_CNNinfo['next_sense_base1'] + df_CNNinfo['next_sense_base2']).map(base_map),
        np.nan
    )

    df_CNNinfo['antisense_motif'] = np.where(
        df_CNNinfo['antisense_base'].str.upper() == 'C',
        (df_CNNinfo['next_antisense_base1'] + df_CNNinfo['next_antisense_base2']).map(base_map),
        np.nan
    )

    # Drop temporary columns
    df_CNNinfo.drop(columns=['next_sense_base1', 'next_sense_base2', 'next_antisense_base1', 'next_antisense_base2'], inplace=True)

    return df_CNNinfo


def calculate_CNN(df):
    """Calculate CNN trimer ratios"""

    df_filtered = df[df['sense_motif'].notna() | df['antisense_motif'].notna()]

    # Create a long-form dataframe for both strands
    trimer_counts = df_filtered.melt(
        value_vars=['sense_motif', 'antisense_motif'],
        var_name='Strand',
        value_name='Trimer_type'
    ).dropna()

    # Map the original rows to the melted dataframe to bring back Nvalid_cov and Nmod
    trimer_counts['Nvalid_cov'] = trimer_counts.index.map(df['Nvalid_cov'])
    trimer_counts['Nmod'] = trimer_counts.index.map(df['Nmod'])

    # Group by Trimer_type and calculate Frequency and methylated_Frequency
    result = trimer_counts.groupby('Trimer_type').agg(
        Frequency=('Nvalid_cov', 'sum'),
        methylated_Frequency=('Nmod', 'sum')
    ).reset_index()

    # calculate ratio
    result['ratio'] = result['methylated_Frequency']/(result['methylated_Frequency']+result['Frequency'])

    return result


def process_CNN_analysis(wholegenome_data):
    """Process CNN trimer analysis and create visualizations"""

    samples_wholegenome = list(wholegenome_data.keys())
    trimer_data = {}

    print("Processing CNN trimer analysis...")

    for sample in samples_wholegenome:
        CNN_info = find_CNN_trimer(wholegenome_data[sample])
        trimer_info = calculate_CNN(CNN_info)
        trimer_data[sample] = trimer_info

    # Create visualizations
    output_dir = "visualization_motif_analysis/CNN_analysis"
    os.makedirs(output_dir, exist_ok=True)

    for sample in samples_wholegenome:
        plt.figure(figsize=(7,4))

        sns.scatterplot(data=trimer_data[sample],
                        x='Frequency',
                        y='ratio',
                        zorder=2)

        for i in range(trimer_data[sample].shape[0]):
            plt.text(x=trimer_data[sample]['Frequency'].iloc[i],
                    y=trimer_data[sample]['ratio'].iloc[i],
                    s=trimer_data[sample]['Trimer_type'].iloc[i])

        plt.ylabel("5mCNN to total CNN ratio")
        plt.xlabel("Total occurrence of CNN")
        plt.title(f"{sample}")
        plt.grid(axis='both', zorder=0, color='lightgrey')

        sns.despine()

        output_path = f"{output_dir}/CNN_trimer_{sample}.png"
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f"Saved CNN analysis: {output_path}")

    return trimer_data


def main():
    """Main function to run the base information and motif analysis"""

    # Setup data
    lpat_genes, samples, contigs = setup_data()
    # lpat_genes, dgat_genes, samples, contigs = setup_data()

    # 1. Load BED files
    print("=== 1. Loading BED Files ===")
    sample_data, BED_colnames = load_bed_files(samples, lpat_genes, contigs)
    # sample_data, BED_colnames = load_bed_files(samples, lpat_genes, dgat_genes, contigs)

    # 2. Generate sense base information
    print("\n=== 2. Base Information for BED Files ===")
    generate_sense_base_info(samples, lpat_genes)
    # generate_sense_base_info(samples, lpat_genes, dgat_genes)

    # 3. Combine base information
    baseinfo_data = combine_base_information(samples, lpat_genes)
    # baseinfo_data = combine_base_information(samples, lpat_genes, dgat_genes)

    # 4. Methylation motif identification
    print("\n=== 3. Methylation Motif Analysis ===")
    motif_data = process_motif_analysis(baseinfo_data, samples, lpat_genes)
    # motif_data = process_motif_analysis(baseinfo_data, samples, lpat_genes, dgat_genes)

    # 5. Create methylation landscape plots
    print("\n=== 4. Creating Methylation Landscape Plots ===")
    create_methylation_landscape_plots(motif_data, samples, lpat_genes, contigs)
    # create_methylation_landscape_plots(motif_data, samples, lpat_genes, dgat_genes, contigs)

    # 6. Whole genome analysis
    print("\n=== 5. Whole Genome Analysis ===")
    generate_whole_genome_base_info()
    wholegenome_data = load_whole_genome_base_info()

    # 7. CN dimer analysis
    print("\n=== 6. CN Dimer Analysis ===")
    dimer_data = process_CN_analysis(wholegenome_data)

    # 8. CNN trimer analysis
    print("\n=== 7. CNN Trimer Analysis ===")
    trimer_data = process_CNN_analysis(wholegenome_data)

    print("\n=== Base Information and Motif Analysis Completed! ===")


if __name__ == "__main__":
    main()
