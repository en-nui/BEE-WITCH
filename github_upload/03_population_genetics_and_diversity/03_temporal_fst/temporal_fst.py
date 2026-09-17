#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Temporal Fst at gene level from SNV frequency trajectories
---------------------------------------------------------

Reads per-colony SNV trajectory files of the form:
  snv_frequency_trajectory_<Colony>.csv

Each file has rows like:
MAG_id,Month,contig_name,pos_in_contig,reference,unique_gene_callers_id,
corresponding_gene_call,A,C,G,T,coverage,polarized_major_allele,
polarized_minor_allele,departure_from_polarized_consensus,
base_pos_in_codon,codon_order_in_gene,codon_number,gene_length

For each Colony × MAG_id × corresponding_gene_call:
  - Compute Fst between adjacent months (Month_{i+1} vs Month_i)
  - Compute Fst between each month and the initial month
  - Compute per-gene averages and normalized Fst

Outputs:
  gene_month_to_month_fst.csv
  gene_vs_initial_month_fst.csv
  gene_average_temporal_fst.csv
"""

import pandas as pd
import numpy as np
from multiprocessing import Pool
import os
import glob

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------

MONTH_ORDER = ['May', 'June', 'July', 'August',
               'September', 'October', 'November',
               'January', 'February']

TRAJ_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories_polymorphic_only"

OUTDIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/temporal_fst_gene_level"
os.makedirs(OUTDIR, exist_ok=True)

N_PROCESSES = 3  # adjust as needed


# ----------------------------------------------------------------------
# Fst helper functions
# ----------------------------------------------------------------------

def calculate_nucleotide_diversity(df):
    """
    Per-sample nucleotide diversity (theta_pi) using A,C,G,T counts and coverage.

    df must have columns: coverage, A, C, G, T

    pi_site = (1 - sum_i p_i^2) * (coverage / (coverage - 1))
    Returns the mean across sites.
    """
    if df.empty:
        return np.nan

    df = df.copy()
    df = df[df['coverage'] > 0]
    if df.empty:
        return np.nan

    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)

    numer = (df['A'] ** 2 + df['C'] ** 2 + df['G'] ** 2 + df['T'] ** 2)
    denom = df['coverage'] ** 2
    df['freq'] = numer / denom

    df['heterozygosity'] = (1.0 - df['freq']) * df['site_cov']

    return df['heterozygosity'].sum() / len(df)


def calculate_heterozygosity_for_combined(df):
    """
    Given a combined data frame with coverage, A, C, G, T for both samples summed,
    compute site-level heterozygosity as above.
    """
    if df.empty:
        df['heterozygosity'] = np.nan
        return df

    df = df.copy()
    df = df[df['coverage'] > 0]
    if df.empty:
        df['heterozygosity'] = np.nan
        return df

    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)
    numer = (df['A'] ** 2 + df['C'] ** 2 + df['G'] ** 2 + df['T'] ** 2)
    denom = df['coverage'] ** 2
    df['freq'] = numer / denom
    df['heterozygosity'] = (1.0 - df['freq']) * df['site_cov']
    return df


def calculate_fst(df1, df2,
                  merge_cols=('contig_name', 'pos_in_contig')):
    """
    Fst between two samples (df1, df2) using heterozygosity-based estimator.

    df1, df2: data frames for a given Colony × MAG_id × gene × Month,
              containing columns:
                - contig_name
                - pos_in_contig
                - coverage, A, C, G, T
    """
    if df1.empty or df2.empty:
        return np.nan, np.nan

    theta1 = calculate_nucleotide_diversity(df1)
    theta2 = calculate_nucleotide_diversity(df2)
    if np.isnan(theta1) or np.isnan(theta2):
        return np.nan, np.nan

    hs = (theta1 + theta2) / 2.0

    df1m = df1[list(merge_cols) + ['coverage', 'A', 'C', 'G', 'T']].copy()
    df2m = df2[list(merge_cols) + ['coverage', 'A', 'C', 'G', 'T']].copy()

    combined_df = pd.merge(
        df1m,
        df2m,
        on=list(merge_cols),
        how='outer',
        suffixes=('_1', '_2')
    )

    for col in ['coverage', 'A', 'C', 'T', 'G']:
        c1 = f"{col}_1"
        c2 = f"{col}_2"
        combined_df[col] = combined_df.get(c1, 0).fillna(0) + \
                           combined_df.get(c2, 0).fillna(0)
        if c1 in combined_df.columns:
            combined_df.drop(columns=[c1], inplace=True)
        if c2 in combined_df.columns:
            combined_df.drop(columns=[c2], inplace=True)

    combined_df = calculate_heterozygosity_for_combined(combined_df)
    if combined_df.empty or combined_df['heterozygosity'].isna().all():
        return np.nan, np.nan

    ht = combined_df['heterozygosity'].sum() / len(combined_df)
    if ht <= 0 or np.isnan(ht):
        return np.nan, np.nan

    fst = (ht - hs) / ht
    fst_method1 = fst
    fst_method2 = fst  # placeholder for a second estimator if desired

    return fst_method1, fst_method2


# ----------------------------------------------------------------------
# Per-MAG / per-gene processing
# ----------------------------------------------------------------------

def process_mag_id_gene_level(mag_id, df):
    """
    Process one MAG_id within a single colony at gene-level.

    df has columns:
      Colony, MAG_id, Month, contig_name, pos_in_contig, A,C,G,T,coverage,
      corresponding_gene_call, ...

    Returns three lists of dicts:
      - gene_month_to_month_results
      - gene_vs_initial_month_results (with normalized Fst)
      - gene_average_results
    """
    print(f"Processing MAG_id: {mag_id} at gene-level on PID {os.getpid()}")
    mag_df = df[df['MAG_id'] == mag_id].copy()
    if mag_df.empty:
        return None, None, None

    colony = mag_df['Colony'].iloc[0]

    # Ensure months are ordered
    mag_df['Month'] = pd.Categorical(
        mag_df['Month'],
        categories=MONTH_ORDER,
        ordered=True
    )

    gene_month_to_month_results = []
    gene_vs_initial_month_results = []
    gene_average_results = []

    unique_genes = mag_df['corresponding_gene_call'].astype(str).unique()

    for gene_id in unique_genes:
        gene_df = mag_df[mag_df['corresponding_gene_call'].astype(str) == gene_id].copy()
        if gene_df.empty:
            continue

        gene_df.sort_values(by='Month', inplace=True)
        unique_months = gene_df['Month'].dropna().unique()

        if len(unique_months) < 2:
            # Not enough time points for this gene
            continue

        initial_month = unique_months[0]
        initial_month_data = gene_df[gene_df['Month'] == initial_month].copy()

        # Adjacent months
        gene_m2m = []
        for i in range(1, len(unique_months)):
            month1 = unique_months[i - 1]
            month2 = unique_months[i]

            data1 = gene_df[gene_df['Month'] == month1].copy()
            data2 = gene_df[gene_df['Month'] == month2].copy()

            fst1, fst2 = calculate_fst(data1, data2)

            gene_m2m.append({
                'Colony': colony,
                'MAG_id': mag_id,
                'corresponding_gene_call': str(gene_id),
                'Comparison': f"{month2} x {month1}",
                'Month1': str(month1),
                'Month2': str(month2),
                'Fst_Method1': fst1,
                'Fst_Method2': fst2
            })

        # Month vs initial
        gene_vs0 = []
        for i in range(1, len(unique_months)):
            month = unique_months[i]
            current_month_data = gene_df[gene_df['Month'] == month].copy()

            fst1, fst2 = calculate_fst(initial_month_data, current_month_data)

            gene_vs0.append({
                'Colony': colony,
                'MAG_id': mag_id,
                'corresponding_gene_call': str(gene_id),
                'Comparison': f"{month} x {initial_month}",
                'Initial_Month': str(initial_month),
                'Month': str(month),
                'Fst_Method1': fst1,
                'Fst_Method2': fst2
            })

        # Averages + normalization for this gene
        df_m2m = pd.DataFrame(gene_m2m) if gene_m2m else None
        df_vs0 = pd.DataFrame(gene_vs0) if gene_vs0 else None

        avg_fst_m2m_1 = df_m2m['Fst_Method1'].mean() if df_m2m is not None else np.nan
        avg_fst_m2m_2 = df_m2m['Fst_Method2'].mean() if df_m2m is not None else np.nan

        gene_vs0_norm = []
        if df_vs0 is not None:
            for _, res in df_vs0.iterrows():
                norm_fst1 = res['Fst_Method1'] / avg_fst_m2m_1 \
                    if avg_fst_m2m_1 not in (0, np.nan) else np.nan
                norm_fst2 = res['Fst_Method2'] / avg_fst_m2m_2 \
                    if avg_fst_m2m_2 not in (0, np.nan) else np.nan

                gene_vs0_norm.append({
                    'Colony': colony,
                    'MAG_id': mag_id,
                    'corresponding_gene_call': str(gene_id),
                    'Comparison': res['Comparison'],
                    'Initial_Month': res['Initial_Month'],
                    'Month': res['Month'],
                    'Fst_Method1': res['Fst_Method1'],
                    'Fst_Method2': res['Fst_Method2'],
                    'Normalized_Fst_Method1': norm_fst1,
                    'Normalized_Fst_Method2': norm_fst2
                })

        avg_result = {
            'Colony': colony,
            'MAG_id': mag_id,
            'corresponding_gene_call': str(gene_id),
            'Average_Fst_Month_to_Month_Method1': avg_fst_m2m_1,
            'Average_Fst_Month_to_Month_Method2': avg_fst_m2m_2,
            'Average_Fst_vs_Initial_Month_Method1':
                df_vs0['Fst_Method1'].mean() if df_vs0 is not None else np.nan,
            'Average_Fst_vs_Initial_Month_Method2':
                df_vs0['Fst_Method2'].mean() if df_vs0 is not None else np.nan,
            'Average_Normalized_Fst_Method1':
                np.nanmean([r['Normalized_Fst_Method1'] for r in gene_vs0_norm])
                if gene_vs0_norm else np.nan,
            'Average_Normalized_Fst_Method2':
                np.nanmean([r['Normalized_Fst_Method2'] for r in gene_vs0_norm])
                if gene_vs0_norm else np.nan
        }

        gene_month_to_month_results.extend(gene_m2m)
        gene_vs_initial_month_results.extend(gene_vs0_norm)
        gene_average_results.append(avg_result)

    return gene_month_to_month_results, gene_vs_initial_month_results, gene_average_results


# ----------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------

def load_and_prepare_trajectory_file(filepath):
    """
    Load one snv_frequency_trajectory_<Colony>.csv and add a Colony column.
    Ensure numeric types for A, C, G, T, coverage.
    """
    df = pd.read_csv(filepath)
    basename = os.path.basename(filepath)
    # snv_frequency_trajectory_100.csv -> "100"
    colony = os.path.splitext(basename)[0].split('_')[-1]
    df['Colony'] = str(colony)

    # numeric counts
    cols_numeric = ['coverage', 'A', 'C', 'G', 'T']
    for col in cols_numeric:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

    df = df[df['coverage'] > 0]

    # ensure string types
    for c in ['MAG_id', 'Month', 'contig_name', 'corresponding_gene_call']:
        if c in df.columns:
            df[c] = df[c].astype(str)

    return df


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    traj_files = sorted(glob.glob(os.path.join(TRAJ_DIR, "snv_frequency_trajectory_*.csv")))
    if not traj_files:
        print(f"No trajectory files found in {TRAJ_DIR}")
        return

    print(f"Found {len(traj_files)} trajectory files:")
    for f in traj_files:
        print("  ", f)

    all_gene_m2m = []
    all_gene_vs0 = []
    all_gene_avg = []

    for traj_file in traj_files:
        print("\n==================================================")
        print("Processing file:", traj_file)
        df = load_and_prepare_trajectory_file(traj_file)

        unique_mag_ids = df['MAG_id'].unique()
        print(f"  Colony {df['Colony'].iloc[0]}: {len(unique_mag_ids)} MAG_id")

        with Pool(processes=N_PROCESSES) as pool:
            results = pool.starmap(
                process_mag_id_gene_level,
                [(mag_id, df) for mag_id in unique_mag_ids]
            )

        for m2m, vs0, avg in results:
            if m2m:
                all_gene_m2m.extend(m2m)
            if vs0:
                all_gene_vs0.extend(vs0)
            if avg:
                all_gene_avg.extend(avg)

    # Save gene-level month-to-month Fst
    if all_gene_m2m:
        df_m2m = pd.DataFrame(all_gene_m2m)
        out_path = os.path.join(OUTDIR, "gene_month_to_month_fst.csv")
        df_m2m.to_csv(out_path, index=False)
        print("\nGene month-to-month Fst results saved to", out_path)
    else:
        print("\nNo gene month-to-month Fst results to save.")

    # Save gene-level vs-initial-month Fst (normalized)
    if all_gene_vs0:
        df_vs0 = pd.DataFrame(all_gene_vs0)
        out_path = os.path.join(OUTDIR, "gene_vs_initial_month_fst.csv")
        df_vs0.to_csv(out_path, index=False)
        print("Gene vs initial month Fst results saved to", out_path)
    else:
        print("\nNo gene vs initial month Fst results to save.")

    # Save gene-level average Fst
    if all_gene_avg:
        df_avg = pd.DataFrame(all_gene_avg)
        out_path = os.path.join(OUTDIR, "gene_average_temporal_fst.csv")
        df_avg.to_csv(out_path, index=False)
        print("Gene average temporal Fst results saved to", out_path)
    else:
        print("\nNo gene average temporal Fst results to save.")

    print("\nFinished processing and saving all gene-level Fst results.")


if __name__ == "__main__":
    main()
