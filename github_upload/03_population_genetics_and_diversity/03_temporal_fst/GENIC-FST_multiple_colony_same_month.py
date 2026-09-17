#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 14:33:05 2025

@author: robinch
"""

import pandas as pd
import numpy as np
from itertools import combinations

# Define the month order (if not already defined)
month_order = ['May', 'June', 'July', 'August', 'September', 'October', 'November', 'January', 'February']

# --- Fst Calculation Functions (same as before) ---
def calculate_nucleotide_diversity(df):
    if df.empty:
        return np.nan
    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)
    df['freq'] = df.apply(
        lambda row: float(
            (row['A'] / row['coverage'])**2 +
            (row['T'] / row['coverage'])**2 +
            (row['C'] / row['coverage'])**2 +
            (row['G'] / row['coverage'])**2
        ) if row['coverage'] > 0 else np.nan,
        axis=1
    )
    df['heterozygosity'] = (1 - df['freq']) * df['site_cov']
    return df['heterozygosity'].sum() / len(df)

def calculate_heterozygosity_for_combined(df):
    if df.empty:
        return df
    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)
    df['freq'] = df.apply(
        lambda row: float(
            (row['A'] / row['coverage'])**2 +
            (row['T'] / row['coverage'])**2 +
            (row['C'] / row['coverage'])**2 +
            (row['G'] / row['coverage'])**2
        ) if row['coverage'] > 0 else np.nan,
        axis=1
    )
    df['heterozygosity'] = (1 - df['freq']) * df['site_cov']
    return df

def calculate_fst(df1, df2, pos_column='pos'):
    if df1.empty or df2.empty:
        return np.nan, np.nan

    theta1 = calculate_nucleotide_diversity(df1)
    theta2 = calculate_nucleotide_diversity(df2)
    hs = (theta1 + theta2) / 2

    combined_df = pd.merge(df1, df2, on=pos_column, how='outer', suffixes=('_1', '_2'))
    cols_to_sum = ['coverage', 'A', 'C', 'T', 'G']
    for col in cols_to_sum:
        combined_df[col] = combined_df[f'{col}_1'].fillna(0) + combined_df[f'{col}_2'].fillna(0)
        if f'{col}_1' in combined_df.columns:
            combined_df.drop(columns=[f'{col}_1', f'{col}_2'], inplace=True)
    combined_df = calculate_heterozygosity_for_combined(combined_df)
    ht = combined_df['heterozygosity'].sum() / len(combined_df) if len(combined_df) > 0 else np.nan

    fst_method1 = (ht - hs) / ht if ht > 0 else np.nan
    fst_method2 = (ht - hs) / ht if ht > 0 else np.nan # User's second method is the same formula

    return fst_method1, fst_method2

# Load the DataFrame for analysis
try:
    multiple_colony_shared_month_df = pd.read_csv('/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genome-wide/multiple_colony_same_month/multiple_colony_shared_month.csv')
    multiple_colony_shared_month_df = multiple_colony_shared_month_df.dropna()
    print("Loaded multiple_colony_shared_month.csv")
    initial_rows = len(multiple_colony_shared_month_df) # Store initial number of rows
except FileNotFoundError:
    print("Error: multiple_colony_shared_month.csv not found. Please ensure the file is in the correct directory.")
    exit()

# Load the retained MAGs data
retained_mags_path = '/home/robinch/projects/BEE-WITCH/21_final_r_vMAG_set/all_colonies_RETAINED_r_vMAGs.csv'
try:
    retained_mags_df = pd.read_csv(retained_mags_path)
    print("Loaded all_colonies_RETAINED_r_vMAGs.csv")
except FileNotFoundError:
    print(f"Error: Retained MAGs file not found at {retained_mags_path}")
    exit()

# --- Filter the main DataFrame based on the retained MAGs data ---
print("Filtering data based on retained MAGs, Colonies, and Months...")
# Merge the two DataFrames to keep only matching rows
multiple_colony_shared_month_df = pd.merge(multiple_colony_shared_month_df,
                                             retained_mags_df[['MAG_id', 'Colony', 'Month']],
                                             on=['MAG_id', 'Colony', 'Month'],
                                             how='inner')
final_rows = len(multiple_colony_shared_month_df) # Store final number of rows
rows_removed = initial_rows - final_rows
print(f"Filtered DataFrame has {final_rows} rows after matching with retained MAGs data.")
print(f"Number of rows removed during filtering: {rows_removed}")

# Dictionary to store all gene-level Fst results
all_gene_fst_results = {}

# Iterate through each unique MAG_id and gene_pos combination
unique_genes = multiple_colony_shared_month_df.groupby(['MAG_id', 'gene_pos']).ngroup().unique()
num_unique_genes = len(unique_genes)
print(f"\nNumber of unique genes (MAG_id, gene_pos) after filtering: {num_unique_genes}")

grouped_by_gene = multiple_colony_shared_month_df.groupby(['MAG_id', 'gene_pos'])

gene_counter = 0
for (mag_id, gene_pos), gene_df in grouped_by_gene:
    gene_counter += 1
    print(f"\nProcessing gene {gene_counter}/{num_unique_genes}: MAG_id = {mag_id}, gene_pos = {gene_pos}")
    gene_results = []

    # Iterate through each unique month present for this gene
    unique_months = sorted(gene_df['Month'].unique(), key=lambda m: month_order.index(m) if m in month_order else len(month_order))
    for month in unique_months:
        month_df = gene_df[gene_df['Month'] == month].copy()
        unique_colonies = month_df['Colony'].unique()

        # Perform pairwise Fst between each Colony combination
        for colony1, colony2 in combinations(unique_colonies, 2):
            df_colony1 = month_df[month_df['Colony'] == colony1].copy()
            df_colony2 = month_df[month_df['Colony'] == colony2].copy()

            # Ensure coverage and allele counts are numeric
            cols_to_numeric = ['coverage', 'A', 'C', 'T', 'G']
            for col in cols_to_numeric:
                if df_colony1[col].dtype == 'object':
                    df_colony1[col] = pd.to_numeric(df_colony1[col], errors='coerce')
                if df_colony2[col].dtype == 'object':
                    df_colony2[col] = pd.to_numeric(df_colony2[col], errors='coerce')

            fst1, fst2 = calculate_fst(df_colony1, df_colony2)

            # Calculate SNV information
            snvs_colony1 = set(df_colony1['pos'].dropna().astype(int))
            snvs_colony2 = set(df_colony2['pos'].dropna().astype(int))
            public_snvs = snvs_colony1.intersection(snvs_colony2)
            private_snvs_colony1 = snvs_colony1 - public_snvs
            private_snvs_colony2 = snvs_colony2 - public_snvs

            gene_results.append({
                'MAG_id': mag_id,
                'gene_pos': gene_pos,
                'Month': month,
                'Colony1': colony1,
                'Colony2': colony2,
                'Fst_Method1': fst1,
                'Fst_Method2': fst2,
                'SNVs_Colony1': len(snvs_colony1),
                'SNVs_Colony2': len(snvs_colony2),
                'Public_SNVs': len(public_snvs),
                'Private_SNVs_Colony1': len(private_snvs_colony1),
                'Private_SNVs_Colony2': len(private_snvs_colony2)
            })

    if gene_results:
        results_df = pd.DataFrame(gene_results)

        # Calculate average Fst per Month for this gene
        average_fst_per_month = results_df.groupby('Month')[['Fst_Method1', 'Fst_Method2']].mean()
        print(f"\n{mag_id}, {gene_pos} - Average Fst per Month:")
        print(average_fst_per_month)

        # Calculate average Fst per Colony comparison for this gene
        results_df['Colony_Comparison'] = results_df.apply(lambda row: tuple(sorted((row['Colony1'], row['Colony2']))), axis=1)
        average_fst_per_colony_pair = results_df.groupby('Colony_Comparison')[['Fst_Method1', 'Fst_Method2']].mean()
        print(f"\n{mag_id}, {gene_pos} - Average Fst per Colony Comparison:")
        print(average_fst_per_colony_pair)

        # Calculate average Fst across all comparisons for this gene
        average_fst_overall = results_df[['Fst_Method1', 'Fst_Method2']].mean()
        print(f"\n{mag_id}, {gene_pos} - Average Fst across all comparisons:")
        print(average_fst_overall)

        all_gene_fst_results[(mag_id, gene_pos)] = {
            'pairwise_fst': results_df.to_dict(orient='records'),
            'average_fst_per_month': average_fst_per_month.to_dict('index'),
            'average_fst_per_colony_pair': average_fst_per_colony_pair.to_dict('index'),
            'average_fst_overall': average_fst_overall.to_dict()
        }

# 1. Save Gene-level Pairwise Fst Results
gene_pairwise_fst_list = []
for (mag_id, gene_pos), results in all_gene_fst_results.items():
    for record in results['pairwise_fst']:
        record['gene_pos'] = gene_pos  # Add gene_pos back to the record
        gene_pairwise_fst_list.append(record)

if gene_pairwise_fst_list:
    gene_pairwise_fst_df = pd.DataFrame(gene_pairwise_fst_list)
    gene_pairwise_fst_df.to_csv('gene_level_pairwise_fst_results.csv', index=False)
    print("\nGene-level Pairwise Fst results saved to gene_level_pairwise_fst_results.csv")
else:
    print("\nNo gene-level pairwise Fst results to save.")

# 2. Save Gene-level Average Fst per Month
gene_average_fst_month_list = []
for (mag_id, gene_pos), results in all_gene_fst_results.items():
    for month, fst_values in results['average_fst_per_month'].items():
        gene_average_fst_month_list.append({
            'MAG_id': mag_id,
            'gene_pos': gene_pos,
            'Month': month,
            'Average_Fst_Method1': fst_values.get('Fst_Method1'),
            'Average_Fst_Method2': fst_values.get('Fst_Method2')
        })

if gene_average_fst_month_list:
    gene_average_fst_month_df = pd.DataFrame(gene_average_fst_month_list)
    gene_average_fst_month_df.to_csv('gene_level_average_fst_per_month.csv', index=False)
    print("Gene-level Average Fst per Month saved to gene_level_average_fst_per_month.csv")
else:
    print("No gene-level average Fst per month results to save.")

# 3. Save Gene-level Average Fst per Colony Pair
gene_average_fst_colony_pair_list = []
for (mag_id, gene_pos), results in all_gene_fst_results.items():
    for colony_pair_tuple, fst_values in results['average_fst_per_colony_pair'].items():
        colony1, colony2 = colony_pair_tuple
        gene_average_fst_colony_pair_list.append({
            'MAG_id': mag_id,
            'gene_pos': gene_pos,
            'Colony1': colony1,
            'Colony2': colony2,
            'Average_Fst_Method1': fst_values.get('Fst_Method1'),
            'Average_Fst_Method2': fst_values.get('Fst_Method2')
        })

if gene_average_fst_colony_pair_list:
    gene_average_fst_colony_pair_df = pd.DataFrame(gene_average_fst_colony_pair_list)
    gene_average_fst_colony_pair_df.to_csv('gene_level_average_fst_per_colony_pair.csv', index=False)
    print("Gene-level Average Fst per Colony Pair saved to gene_level_average_fst_per_colony_pair.csv")
else:
    print("No gene-level average Fst per colony pair results to save.")

# 4. Save Gene-level Overall Average Fst
gene_overall_average_fst_list = []
for (mag_id, gene_pos), results in all_gene_fst_results.items():
    overall_fst_values = results['average_fst_overall']
    gene_overall_average_fst_list.append({
        'MAG_id': mag_id,
        'gene_pos': gene_pos,
        'Average_Fst_Method1': overall_fst_values.get('Fst_Method1'),
        'Average_Fst_Method2': overall_fst_values.get('Fst_Method2')
    })

if gene_overall_average_fst_list:
    gene_overall_average_fst_df = pd.DataFrame(gene_overall_average_fst_list)
    gene_overall_average_fst_df.to_csv('gene_level_overall_average_fst.csv', index=False)
    print("Gene-level Overall Average Fst saved to gene_level_overall_average_fst.csv")
else:
    print("No gene-level overall average Fst results to save.")

print("\nFinished processing and saving all gene-level results for multiple colonies (same month).")