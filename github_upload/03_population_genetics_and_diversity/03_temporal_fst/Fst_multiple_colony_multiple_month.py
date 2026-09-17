#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 2025

@author: robinch
"""

import os
import pandas as pd
import numpy as np
from itertools import combinations

# Define the month order
month_order = ['May', 'June', 'July', 'August', 'September', 'October', 'November', 'January', 'February']

# --- Fst Calculation Functions (modified to return DataFrames) ---
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
        return (np.nan, np.nan), df1, df2

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

    return (fst_method1, fst_method2), df1, df2

# Load the DataFrame
try:
    multiple_colony_multiple_months_df = pd.read_csv('/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/multiple_colony_multiple_months.csv')
    print("Loaded multiple_colony_multiple_months.csv")
except FileNotFoundError:
    print("Error: multiple_colony_multiple_months.csv not found. Please ensure the file is in the correct directory.")
    exit()

# Identify eligible MAG_ids
eligible_mag_ids = []
grouped_by_mag = multiple_colony_multiple_months_df.groupby('MAG_id')

for mag_id, group in grouped_by_mag:
    unique_colonies = group['Colony'].unique()
    if len(unique_colonies) >= 2:
        is_eligible = True
        for col1, col2 in combinations(unique_colonies, 2):
            months_col1 = set(group[group['Colony'] == col1]['Month'].unique())
            months_col2 = set(group[group['Colony'] == col2]['Month'].unique())
            shared_months = months_col1.intersection(months_col2)
            if len(shared_months) < 2:
                is_eligible = False
                break
        if is_eligible:
            eligible_mag_ids.append(mag_id)

print(f"\nNumber of eligible MAG_ids for analysis: {len(eligible_mag_ids)}")

all_temporal_fst_results = []
all_between_colony_fst_results = []
all_normalized_fst_results = []
average_between_colony_fst_per_mag = {}

# First, calculate all between-colony Fst for each eligible MAG_id
for mag_id in eligible_mag_ids:
    print(f"\nCalculating between-colony Fst for MAG_id: {mag_id}")
    mag_id_df = multiple_colony_multiple_months_df[multiple_colony_multiple_months_df['MAG_id'] == mag_id].copy()
    unique_colonies = sorted(mag_id_df['Colony'].unique())
    unique_months = sorted([month for month in month_order if month in mag_id_df['Month'].unique()], key=lambda m: month_order.index(m))

    between_colony_fsts_mag = []
    for month in unique_months:
        month_df = mag_id_df[mag_id_df['Month'] == month].copy()
        colonies_in_month = sorted(month_df['Colony'].unique())

        if len(colonies_in_month) < 2:
            continue

        for colony1, colony2 in combinations(colonies_in_month, 2):
            df_colony1 = month_df[month_df['Colony'] == colony1].copy()
            df_colony2 = month_df[month_df['Colony'] == colony2].copy()

            cols_to_numeric = ['coverage', 'A', 'C', 'T', 'G']
            for col in cols_to_numeric:
                if df_colony1[col].dtype == 'object':
                    df_colony1[col] = pd.to_numeric(df_colony1[col], errors='coerce')
                if df_colony2[col].dtype == 'object':
                    df_colony2[col] = pd.to_numeric(df_colony2[col], errors='coerce')

            (fst1, fst2), res_df1, res_df2 = calculate_fst(df_colony1, df_colony2)

            unique_snvs_c1 = set(res_df1['pos'].dropna().astype(int))
            unique_snvs_c2 = set(res_df2['pos'].dropna().astype(int))
            shared_snvs_c = unique_snvs_c1.intersection(unique_snvs_c2)

            all_between_colony_fst_results.append({
                'MAG_id': mag_id,
                'Month': month,
                'Colony1': colony1,
                'Colony2': colony2,
                'Fst_Method1': fst1,
                'Fst_Method2': fst2,
                'SNVs_Colony1': len(unique_snvs_c1),
                'SNVs_Colony2': len(unique_snvs_c2),
                'Shared_SNVs': len(shared_snvs_c),
                'Unique_to_Colony1': len(unique_snvs_c1 - shared_snvs_c),
                'Unique_to_Colony2': len(unique_snvs_c2 - shared_snvs_c)
            })
            if not np.isnan(fst1):
                between_colony_fsts_mag.append(fst1)

    average_between_colony_fst_per_mag[mag_id] = np.nanmean(between_colony_fsts_mag) if between_colony_fsts_mag else np.nan

# Now, calculate temporal Fst and normalize using the pre-calculated averages
for mag_id in eligible_mag_ids:
    print(f"\nCalculating temporal and normalized Fst for MAG_id: {mag_id}")
    mag_id_df = multiple_colony_multiple_months_df[multiple_colony_multiple_months_df['MAG_id'] == mag_id].copy()
    unique_colonies = sorted(mag_id_df['Colony'].unique())
    unique_months = sorted([month for month in month_order if month in mag_id_df['Month'].unique()], key=lambda m: month_order.index(m))

    for colony in unique_colonies:
        colony_df = mag_id_df[mag_id_df['Colony'] == colony].copy()
        present_months = sorted([month for month in month_order if month in colony_df['Month'].unique()], key=lambda m: month_order.index(m))

        if not present_months or len(present_months) < 2:
            continue

        initial_month = present_months[0]
        df_initial_month = colony_df[colony_df['Month'] == initial_month].copy()

        for month in present_months:
            if month == initial_month:
                continue

            df_month = colony_df[colony_df['Month'] == month].copy()

            cols_to_numeric = ['coverage', 'A', 'C', 'T', 'G']
            for col in cols_to_numeric:
                if df_initial_month[col].dtype == 'object':
                    df_initial_month[col] = pd.to_numeric(df_initial_month[col], errors='coerce')
                if df_month[col].dtype == 'object':
                    df_month[col] = pd.to_numeric(df_month[col], errors='coerce')

            (fst1, fst2), res_df_m1, res_df_m2 = calculate_fst(df_initial_month, df_month)

            unique_snvs_m1 = set(res_df_m1['pos'].dropna().astype(int))
            unique_snvs_m2 = set(res_df_m2['pos'].dropna().astype(int))
            shared_snvs_m = unique_snvs_m1.intersection(unique_snvs_m2)

            temporal_result = {
                'MAG_id': mag_id,
                'Colony': colony,
                'Comparison': f'{initial_month} vs {month}',
                'Month1': initial_month,
                'Month2': month,
                'Fst_Method1': fst1,
                'Fst_Method2': fst2,
                'SNVs_Month1': len(unique_snvs_m1),
                'SNVs_Month2': len(unique_snvs_m2),
                'Shared_SNVs': len(shared_snvs_m),
                'Unique_to_Month1': len(unique_snvs_m1 - shared_snvs_m),
                'Unique_to_Month2': len(unique_snvs_m2 - shared_snvs_m)
            }
            all_temporal_fst_results.append(temporal_result)

            avg_between_colony_fst = average_between_colony_fst_per_mag.get(mag_id)
            normalized_fst1 = np.nan
            normalized_fst2 = np.nan
            normalized_flag = False

            if avg_between_colony_fst is not None and avg_between_colony_fst > 0:
                normalized_flag = True
                if fst1 is not None and fst1 >= 0:
                    normalized_fst1 = fst1 / avg_between_colony_fst
                else:
                    normalized_fst1 = 0
                if fst2 is not None and fst2 >= 0:
                    normalized_fst2 = fst2 / avg_between_colony_fst
                else:
                    normalized_fst2 = 0
            elif avg_between_colony_fst is not None and avg_between_colony_fst <= 0:
                normalized_fst1 = fst1
                normalized_fst2 = fst2

            normalized_result = temporal_result.copy()
            normalized_result.update({
                'Average_Between_Colony_Fst_Method1': avg_between_colony_fst,
                'Normalized_Fst_Method1': normalized_fst1,
                'Normalized_Fst_Method2': normalized_fst2,
                'normalized_flag': normalized_flag
            })
            all_normalized_fst_results.append(normalized_result)

# Save results to CSVs
if all_temporal_fst_results:
    df_temporal_fst = pd.DataFrame(all_temporal_fst_results)
    df_temporal_fst.to_csv('temporal_fst_multiple_colony.csv', index=False)
    print("\nTemporal Fst results saved to temporal_fst_multiple_colony.csv")
else:
    print("\nNo temporal Fst results to save.")

if all_between_colony_fst_results:
    df_between_colony_fst = pd.DataFrame(all_between_colony_fst_results)
    df_between_colony_fst.to_csv('between_colony_fst_multiple_colony.csv', index=False)
    print("Between-Colony Fst results saved to between_colony_fst_multiple_colony.csv")
else:
    print("\nNo between-colony Fst results to save.")

if all_normalized_fst_results:
    df_normalized_fst = pd.DataFrame(all_normalized_fst_results)
    df_normalized_fst.to_csv('normalized_fst_multiple_colony.csv', index=False)
    print("Normalized Fst results saved to normalized_fst_multiple_colony.csv")
else:
    print("\nNo normalized Fst results to save.")

print("\nFinished processing MAG_id in multiple colony and multiple month with updated reporting.")