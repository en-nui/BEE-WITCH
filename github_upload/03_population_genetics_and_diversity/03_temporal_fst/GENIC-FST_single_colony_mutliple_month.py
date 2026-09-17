#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 16:52:06 2025

@author: robinch


@author: robinch
"""

import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import os
from tqdm import tqdm # Using tqdm for progress bars during chunk processing

# Define the month order
month_order = ['May', 'June', 'July', 'August', 'September', 'October', 'November', 'January', 'February']

# Define file paths
single_colony_data_path = '/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genome-wide/single_colony_multiple_month/single_colony_multiple_months.csv'
retained_mags_path = '/home/robinch/projects/BEE-WITCH/21_final_r_vMAG_set/all_colonies_RETAINED_r_vMAGs.csv'

# --- Fst Calculation Functions (same as before, modified to return nucleotide diversity) ---
def calculate_nucleotide_diversity(df):
    if df.empty:
        return np.nan
    # Ensure coverage is numeric before calculation
    df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce')
    # Filter out rows with coverage <= 1 to avoid division by zero
    df_filtered = df[df['coverage'] > 1].copy()

    if df_filtered.empty:
        return np.nan

    df_filtered['site_cov'] = df_filtered['coverage'] / (df_filtered['coverage'] - 1)
    df_filtered['freq'] = df_filtered.apply(
        lambda row: float(
            (pd.to_numeric(row['A'], errors='coerce') / row['coverage'])**2 +
            (pd.to_numeric(row['T'], errors='coerce') / row['coverage'])**2 +
            (pd.to_numeric(row['C'], errors='coerce') / row['coverage'])**2 +
            (pd.to_numeric(row['G'], errors='coerce') / row['coverage'])**2
        ),
        axis=1
    )
    df_filtered['heterozygosity'] = (1 - df_filtered['freq']) * df_filtered['site_cov']
    # Sum only non-nan heterozygosity values and divide by the number of such sites
    return df_filtered['heterozygosity'].sum() / (~np.isnan(df_filtered['heterozygosity'])).sum()


def calculate_heterozygosity_for_combined(df):
    if df.empty:
        return df
     # Ensure coverage and allele counts are numeric before calculation
    for col in ['coverage', 'A', 'C', 'T', 'G']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Filter out rows with coverage <= 1 to avoid division by zero
    df_filtered = df[df['coverage'] > 1].copy()

    if df_filtered.empty:
        return pd.DataFrame(columns=df.columns.tolist() + ['site_cov', 'freq', 'heterozygosity']) # Return empty df with new columns


    df_filtered['site_cov'] = df_filtered['coverage'] / (df_filtered['coverage'] - 1)
    df_filtered['freq'] = df_filtered.apply(
        lambda row: float(
            (row['A'] / row['coverage'])**2 +
            (row['T'] / row['coverage'])**2 +
            (row['C'] / row['coverage'])**2 +
            (row['G'] / row['coverage'])**2
        ),
        axis=1
    )
    df_filtered['heterozygosity'] = (1 - df_filtered['freq']) * df_filtered['site_cov']
    return df_filtered

def calculate_fst(df1, df2, pos_column='pos'):
    if df1.empty or df2.empty:
        return np.nan, np.nan, np.nan, np.nan

    # Ensure pos is integer for merging
    for df in [df1, df2]:
        if pos_column in df.columns:
             df[pos_column] = pd.to_numeric(df[pos_column], errors='coerce') # Convert to numeric first
             df.dropna(subset=[pos_column], inplace=True) # Drop rows where pos is NaN after conversion
             df[pos_column] = df[pos_column].astype(int) # Convert to int


    theta1 = calculate_nucleotide_diversity(df1)
    theta2 = calculate_nucleotide_diversity(df2)
    hs = (theta1 + theta2) / 2

    combined_df = pd.merge(df1, df2, on=pos_column, how='outer', suffixes=('_1', '_2'))
    cols_to_sum = ['coverage', 'A', 'C', 'T', 'G']
    for col in cols_to_sum:
        combined_df[col] = combined_df[f'{col}_1'].fillna(0) + combined_df[f'{col}_2'].fillna(0)
        if f'{col}_1' in combined_df.columns:
            combined_df.drop(columns=[f'{col}_1', f'{col}_2'], inplace=True)

    # Recalculate heterozygosity on the combined dataframe
    combined_df = calculate_heterozygosity_for_combined(combined_df)

    # Calculate ht from the combined dataframe's heterozygosity
    ht = combined_df['heterozygosity'].sum() / (~np.isnan(combined_df['heterozygosity'])).sum() if (~np.isnan(combined_df['heterozygosity'])).sum() > 0 else np.nan


    fst_method1 = (ht - hs) / ht if ht > 0 else np.nan
    fst_method2 = (ht - hs) / ht if ht > 0 else np.nan # User's second method is the same formula

    return fst_method1, fst_method2, theta1, theta2


# Function to process data for a single gene (receives the entire DataFrame)
def process_gene(gene_tuple, df):
    mag_id, gene_pos = gene_tuple
    print(f"Processing gene: MAG_id = {mag_id}, gene_pos = {gene_pos} on process {os.getpid()}")

    # Filter the already loaded DataFrame for the specific gene
    gene_df = df[(df['MAG_id'] == mag_id) & (df['gene_pos'] == gene_pos)].copy()


    if gene_df.empty:
        print(f"  No data found for MAG_id = {mag_id}, gene_pos = {gene_pos} after initial filtering.")
        return None, None, None


    # Colony information is already filtered by the initial merge
    colony = gene_df['Colony'].unique()[0]


    gene_df['Month'] = pd.Categorical(gene_df['Month'], categories=month_order, ordered=True)
    gene_df_sorted = gene_df.sort_values(by='Month')
    unique_months = gene_df_sorted['Month'].unique()


    if len(unique_months) < 2:
        print(f"  Insufficient unique months ({len(unique_months)}) after filtering for Fst calculation for MAG_id = {mag_id}, gene_pos = {gene_pos}.")
        return None, None, None

    month_to_month_results = []
    vs_initial_month_results = []

    initial_month = unique_months[0]
    initial_month_data = gene_df_sorted[gene_df_sorted['Month'] == initial_month].copy()

    # Calculate Fst for Month to preceding month
    for i in range(1, len(unique_months)):
        month1 = unique_months[i-1]
        month2 = unique_months[i]
        data1 = gene_df_sorted[gene_df_sorted['Month'] == month1].copy()
        data2 = gene_df_sorted[gene_df_sorted['Month'] == month2].copy()


        fst1, fst2, theta1, theta2 = calculate_fst(data1, data2)

        # Calculate SNV information
        snvs_month1 = set(data1['pos'].dropna().astype(int))
        snvs_month2 = set(data2['pos'].dropna().astype(int))
        shared_snvs = snvs_month1.intersection(snvs_month2)

        month_to_month_results.append({
            'MAG_id': mag_id,
            'gene_pos': gene_pos,
            'Comparison': f"{month2} x {month1}",
            'Month1': month1,
            'Month2': month2,
            'Fst_Method1': fst1,
            'Fst_Method2': fst2,
            'Nucleotide_Diversity_Month1': theta1,
            'Nucleotide_Diversity_Month2': theta2,
            'SNVs_Month1': len(snvs_month1),
            'SNVs_Month2': len(snvs_month2),
            'Shared_SNVs': len(shared_snvs),
            'Unique_to_Month1': len(snvs_month1 - shared_snvs),
            'Unique_to_Month2': len(snvs_month2 - shared_snvs)
        })

    # Calculate Fst for Month to initial month
    for i in range(1, len(unique_months)):
        month = unique_months[i]
        current_month_data = gene_df_sorted[gene_df_sorted['Month'] == month].copy()


        fst1, fst2, theta_initial, theta_current = calculate_fst(initial_month_data, current_month_data)


        # Calculate SNV information
        snvs_initial = set(initial_month_data['pos'].dropna().astype(int))
        snvs_current = set(current_month_data['pos'].dropna().astype(int))
        shared_snvs = snvs_initial.intersection(snvs_current)


        vs_initial_month_results.append({
            'MAG_id': mag_id,
            'gene_pos': gene_pos,
            'Comparison': f"{month} x {initial_month}",
            'Month1': initial_month,
            'Month2': month,
            'Fst_Method1': fst1,
            'Fst_Method2': fst2,
            'Nucleotide_Diversity_Month1': theta_initial,
            'Nucleotide_Diversity_Month2': theta_current,
            'SNVs_Month1': len(snvs_initial),
            'SNVs_Month2': len(snvs_current),
            'Shared_SNVs': len(shared_snvs),
            'Unique_to_Month1': len(snvs_initial - shared_snvs),
            'Unique_to_Month2': len(snvs_current - shared_snvs)
        })

    avg_fst_month_to_month_1 = pd.DataFrame(month_to_month_results)['Fst_Method1'].mean() if month_to_month_results else np.nan
    avg_fst_month_to_month_2 = pd.DataFrame(month_to_month_results)['Fst_Method2'].mean() if month_to_month_results else np.nan

    normalized_fst_results = []
    # Note: Average Between Colony Fst is not calculated in this script (single colony)
    # We will skip the normalization calculation as it requires avg_between_colony_fst

    average_results = {
        'MAG_id': mag_id,
        'gene_pos': gene_pos,
        'Average_Fst_Month_to_Month_Method1': avg_fst_month_to_month_1,
        'Average_Fst_Month_to_Month_Method2': avg_fst_month_to_month_2,
        'Average_Fst_vs_Initial_Month_Method1': pd.DataFrame(vs_initial_month_results)['Fst_Method1'].mean() if vs_initial_month_results else np.nan,
        'Average_Fst_vs_Initial_Month_Method2': pd.DataFrame(vs_initial_month_results)['Fst_Method2'].mean() if vs_initial_month_results else np.nan,
         # Normalized Fst averages are not applicable here
        'Average_Normalized_Fst_Method1': np.nan,
        'Average_Normalized_Fst_Method2': np.nan
    }


    return month_to_month_results, vs_initial_month_results, average_results


# --- Main Execution Block ---
if __name__ == "__main__":
    print("Starting the Fst calculation script (Single Colony) - Memory-efficient initial filtering.")

    # Load the retained MAGs data once
    print("Loading all_colonies_RETAINED_r_vMAGs.csv...")
    try:
        retained_mags_df = pd.read_csv(retained_mags_path)
        print("Loaded all_colonies_RETAINED_r_vMAGs.csv successfully.")
    except FileNotFoundError:
        print(f"Error: Retained MAGs file not found at {retained_mags_path}")
        exit()

    # Create a set of (MAG_id, Colony, Month) tuples from retained data for efficient filtering
    print("Creating set of retained (MAG_id, Colony, Month) tuples...")
    # Ensure consistent data types before creating the set
    retained_mags_df['MAG_id'] = retained_mags_df['MAG_id'].astype(str).str.strip()
    retained_mags_df['Colony'] = retained_mags_df['Colony'].astype(str).str.strip()
    retained_mags_df['Month'] = retained_mags_df['Month'].astype(str).str.strip()

    retained_tuples = set(retained_mags_df[['MAG_id', 'Colony', 'Month']].itertuples(index=False, name=None))
    print(f"Created set with {len(retained_tuples)} unique retained combinations.")
    del retained_mags_df # Free up memory after creating the set


    # --- Memory-efficient Initial Filtering using chunking and set membership ---
    print("Performing memory-efficient initial filtering...")
    filtered_chunks = []
    chunk_size_data = 100000 # Adjust this chunk size based on your RAM

    try:
        # Get the total number of rows for the progress bar (optional but helpful)
        total_rows_initial = sum(1 for row in open(single_colony_data_path, 'r')) - 1 # Subtract 1 for header

        for i, chunk in enumerate(tqdm(pd.read_csv(single_colony_data_path, sep=',', header=0, chunksize=chunk_size_data),
                                       total=(total_rows_initial // chunk_size_data) + 1, # Estimate total chunks
                                       desc="Filtering Chunks")):

            # Ensure columns for filtering are of consistent types (string) and strip whitespace
            for col in ['MAG_id', 'Colony', 'Month']:
                if col in chunk.columns:
                    chunk[col] = chunk[col].astype(str).str.strip()
                else:
                    print(f"Warning: Filtering column '{col}' not found in a data chunk.")
                    # Handle this warning - maybe skip the chunk or log an error

            # Filter the chunk using the set of retained tuples
            chunk_filtered = chunk[
                chunk[['MAG_id', 'Colony', 'Month']].apply(tuple, axis=1).isin(retained_tuples)
            ].copy() # Use .copy() to avoid SettingWithCopyWarning

            filtered_chunks.append(chunk_filtered)

        # Concatenate all filtered chunks into a single DataFrame
        single_colony_multiple_month_df = pd.concat(filtered_chunks, ignore_index=True)

        initial_rows = total_rows_initial # Use the estimated total rows
        final_rows = len(single_colony_multiple_month_df)
        rows_removed = initial_rows - final_rows
        print(f"Filtered DataFrame has {final_rows} rows after initial filtering.")
        print(f"Number of rows removed during initial filtering: {rows_removed}")
        del filtered_chunks # Free up memory


    except FileNotFoundError:
        print(f"Error: {single_colony_data_path} not found during chunked filtering.")
        exit()
    except Exception as e:
        print(f"An error occurred during initial chunked filtering: {e}")
        exit()


    # Group by MAG_id and then get unique gene positions from the filtered DataFrame
    print("Identifying unique genes from the filtered DataFrame...")
    if single_colony_multiple_month_df.empty:
        print("Filtered DataFrame is empty. No genes to process.")
        unique_genes = []
    else:
        grouped_by_mag = single_colony_multiple_month_df.groupby('MAG_id')['gene_pos'].unique()
        unique_genes = []
        for mag_id, gene_positions in grouped_by_mag.items():
            for gene_pos in gene_positions:
                unique_genes.append((mag_id, gene_pos))

        num_total_genes = len(unique_genes)
        print(f"Found {num_total_genes} unique genes to process.")


    num_processes = 1 # You can adjust this based on your CPU cores
    print(f"Using {num_processes} CPUs for processing genes.")

    all_month_to_month_data = []
    all_vs_initial_month_data = []
    all_average_data = []

    if unique_genes: # Only start multiprocessing if there are genes to process
        print("Starting gene processing using multiprocessing...")
        # Use pool.starmap and pass the entire filtered DataFrame to each process
        with Pool(processes=num_processes) as pool:
            # Use tqdm to show progress during multiprocessing
            results = list(tqdm(pool.starmap(process_gene, [(gene, single_colony_multiple_month_df) for gene in unique_genes]),
                                total=num_total_genes, desc="Processing Genes"))


            for month_to_month, vs_initial, averages in results:
                if month_to_month:
                    all_month_to_month_data.extend(month_to_month)
                if vs_initial:
                    all_vs_initial_month_data.extend(vs_initial)
                if averages:
                    all_average_data.append(averages)

        print("Finished processing all genes.")
    else:
        print("No unique genes found after filtering. Skipping gene processing.")


    # Save all month-to-month results
    print("Saving month-to-month Fst results...")
    if all_month_to_month_data:
        df_month_to_month = pd.DataFrame(all_month_to_month_data)
        df_month_to_month.to_csv("gene_level_month_to_month_fst_single_colony.csv", index=False)
        print("\nAll gene-level month-to-month Fst results saved to gene_level_month_to_month_fst_single_colony.csv")
    else:
        print("\nNo gene-level month-to-month Fst results to save.")

    # Save all vs initial month results (including normalized Fst)
    print("Saving vs initial month Fst results...")
    if all_vs_initial_month_data:
        df_vs_initial_month = pd.DataFrame(all_vs_initial_month_data)
        df_vs_initial_month.to_csv("gene_level_vs_initial_month_fst_single_colony.csv", index=False)
        print("All gene-level vs initial month Fst results saved to gene_level_vs_initial_month_fst_single_colony.csv")
    else:
        print("\nNo gene-level vs initial month Fst results to save.")

    # Save all average results
    print("Saving average temporal Fst results...")
    if all_average_data:
        df_average = pd.DataFrame(all_average_data)
        df_average.to_csv("gene_level_average_temporal_fst_single_colony.csv", index=False)
        print("All gene-level average temporal Fst results saved to gene_level_average_temporal_fst_single_colony.csv")
    else:
        print("\nNo gene-level average temporal Fst results to save.")

    print("\nFinished processing and saving all gene-level results for single colony.")