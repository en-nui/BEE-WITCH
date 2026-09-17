#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 2025

@author: robinch
"""

import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import os

# Define the month order
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

def process_mag_id(mag_id, df):
    print(f"Processing MAG_id: {mag_id} on process {os.getpid()}")
    mag_id_df = df[df['MAG_id'] == mag_id].copy()
    colony = mag_id_df['Colony'].unique()[0]

    mag_id_df['Month'] = pd.Categorical(mag_id_df['Month'], categories=month_order, ordered=True)
    mag_id_df_sorted = mag_id_df.sort_values(by='Month')
    unique_months = mag_id_df_sorted['Month'].unique()

    if len(unique_months) < 2:
        print(f"  Insufficient months for Fst calculation for {mag_id}.")
        return None, None, None

    month_to_month_results = []
    vs_initial_month_results = []

    initial_month_data = mag_id_df_sorted[mag_id_df_sorted['Month'] == unique_months[0]].copy()

    # Calculate Fst for Month to preceding month
    for i in range(1, len(unique_months)):
        month1 = unique_months[i-1]
        month2 = unique_months[i]
        data1 = mag_id_df_sorted[mag_id_df_sorted['Month'] == month1].copy()
        data2 = mag_id_df_sorted[mag_id_df_sorted['Month'] == month2].copy()

        cols_to_numeric = ['coverage', 'A', 'C', 'T', 'G']
        for col in cols_to_numeric:
            if data1[col].dtype == 'object':
                data1[col] = pd.to_numeric(data1[col], errors='coerce')
            if data2[col].dtype == 'object':
                data2[col] = pd.to_numeric(data2[col], errors='coerce')

        fst1, fst2 = calculate_fst(data1, data2)
        month_to_month_results.append({
            'MAG_id': mag_id,
            'Comparison': f"{month2} x {month1}",
            'Fst_Method1': fst1,
            'Fst_Method2': fst2
        })

    # Calculate Fst for Month to initial month
    for i in range(1, len(unique_months)):
        month = unique_months[i]
        current_month_data = mag_id_df_sorted[mag_id_df_sorted['Month'] == month].copy()

        cols_to_numeric = ['coverage', 'A', 'C', 'T', 'G']
        for col in cols_to_numeric:
            if initial_month_data[col].dtype == 'object':
                initial_month_data[col] = pd.to_numeric(initial_month_data[col], errors='coerce')
            if current_month_data[col].dtype == 'object':
                current_month_data[col] = pd.to_numeric(current_month_data[col], errors='coerce')

        fst1, fst2 = calculate_fst(initial_month_data, current_month_data)
        vs_initial_month_results.append({
            'MAG_id': mag_id,
            'Comparison': f"{month} x {unique_months[0]}",
            'Fst_Method1': fst1,
            'Fst_Method2': fst2
        })

    avg_fst_month_to_month_1 = pd.DataFrame(month_to_month_results)['Fst_Method1'].mean() if month_to_month_results else np.nan
    avg_fst_month_to_month_2 = pd.DataFrame(month_to_month_results)['Fst_Method2'].mean() if month_to_month_results else np.nan

    normalized_fst_results = []
    for res in vs_initial_month_results:
        norm_fst1 = res['Fst_Method1'] / avg_fst_month_to_month_1 if not np.isnan(avg_fst_month_to_month_1) and avg_fst_month_to_month_1 != 0 else np.nan
        norm_fst2 = res['Fst_Method2'] / avg_fst_month_to_month_2 if not np.isnan(avg_fst_month_to_month_2) and avg_fst_month_to_month_2 != 0 else np.nan
        normalized_fst_results.append({
            'MAG_id': mag_id,
            'Comparison': res['Comparison'],
            'Fst_Method1': res['Fst_Method1'],
            'Fst_Method2': res['Fst_Method2'],
            'Normalized_Fst_Method1': norm_fst1,
            'Normalized_Fst_Method2': norm_fst2
        })

    average_results = {
        'MAG_id': mag_id,
        'Average_Fst_Month_to_Month_Method1': avg_fst_month_to_month_1,
        'Average_Fst_Month_to_Month_Method2': avg_fst_month_to_month_2,
        'Average_Fst_vs_Initial_Month_Method1': pd.DataFrame(vs_initial_month_results)['Fst_Method1'].mean() if vs_initial_month_results else np.nan,
        'Average_Fst_vs_Initial_Month_Method2': pd.DataFrame(vs_initial_month_results)['Fst_Method2'].mean() if vs_initial_month_results else np.nan,
        'Average_Normalized_Fst_Method1': pd.DataFrame(normalized_fst_results)['Normalized_Fst_Method1'].mean() if normalized_fst_results else np.nan,
        'Average_Normalized_Fst_Method2': pd.DataFrame(normalized_fst_results)['Normalized_Fst_Method2'].mean() if normalized_fst_results else np.nan
    }

    return month_to_month_results, normalized_fst_results, average_results

# Load the DataFrame
try:
    single_colony_multiple_month_df = pd.read_csv('/home/robinch/projects/BEE-WITCH/22_anvio/single_colony_multiple_months.csv')
    print("Loaded single_colony_multiple_months.csv")
except FileNotFoundError:
    print("Error: single_colony_multiple_months.csv not found. Please ensure the file is in the correct directory.")
    exit()

if __name__ == "__main__":
    unique_mag_ids = single_colony_multiple_month_df['MAG_id'].unique()
    num_processes = 2
    print(f"Using {num_processes} CPUs for processing.")

    all_month_to_month_data = []
    all_vs_initial_month_data = []
    all_average_data = []

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(process_mag_id, [(mag_id, single_colony_multiple_month_df) for mag_id in unique_mag_ids])

        for month_to_month, vs_initial, averages in results:
            if month_to_month:
                all_month_to_month_data.extend(month_to_month)
            if vs_initial:
                all_vs_initial_month_data.extend(vs_initial)
            if averages:
                all_average_data.append(averages)

    # Save all month-to-month results
    if all_month_to_month_data:
        df_month_to_month = pd.DataFrame(all_month_to_month_data)
        df_month_to_month.to_csv("all_month_to_month_fst.csv", index=False)
        print("\nAll month-to-month Fst results saved to all_month_to_month_fst.csv")
    else:
        print("\nNo month-to-month Fst results to save.")

    # Save all vs initial month results (including normalized Fst)
    if all_vs_initial_month_data:
        df_vs_initial_month = pd.DataFrame(all_vs_initial_month_data)
        df_vs_initial_month.to_csv("all_vs_initial_month_fst.csv", index=False)
        print("All vs initial month Fst results saved to all_vs_initial_month_fst.csv")
    else:
        print("\nNo vs initial month Fst results to save.")

    # Save all average results
    if all_average_data:
        df_average = pd.DataFrame(all_average_data)
        df_average.to_csv("all_average_temporal_fst.csv", index=False)
        print("All average temporal Fst results saved to all_average_temporal_fst.csv")
    else:
        print("\nNo average temporal Fst results to save.")

    print("\nFinished processing and saving all results.")