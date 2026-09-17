#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import glob
from tqdm import tqdm
from scipy.stats import binom, poisson
from statsmodels.sandbox.stats.multicomp import multipletests
import traceback

# --- Configuration ---
# MODIFIED: Input is now a directory
SNV_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories/"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/"

# Filtering thresholds (same as before)
PI_THRESHOLD = 1e-3
MIN_FREQ_THRESHOLD = 0.2 # Max frequency allowed for filtering low-freq sites
LOW_FREQ_SWEEP = 0.2 # Lower bound for sweep definition
HIGH_FREQ_SWEEP = 0.7 # Upper bound for sweep definition
COVERAGE_RATIO_FOLD_CHANGE = 2.0

# Significance thresholds (same as before)
FDR_THRESHOLD = 0.1
BH_ALPHA = 0.05

# Month order (crucial for consecutive pairs)
MONTH_ORDER = ["May", "June", "July", "August", "September", "October", "November", "January", "February"]

# --- Helper Function (Unchanged) ---
def calculate_pi(cov1, cov2, alt1, alt2, freq_thresh_low, freq_thresh_high):
    """Calculates Pi(t1, t2), the probability of observing a large shift under H0."""
    total_cov = cov1 + cov2
    if total_cov <= 0:
        return 1.0 # Cannot calculate f_bar, assume max probability (will be filtered)

    total_alt = alt1 + alt2
    f_bar = total_alt / total_cov
    f_bar = np.clip(f_bar, 1e-9, 1 - 1e-9) # Clip f_bar to avoid p=0 or p=1 for binomial

    k1_low = np.floor(freq_thresh_low * cov1)
    k2_high_ref_interpretation = np.floor( (1.0 - freq_thresh_high) * cov2) # = floor(0.3 * cov2)
    prob_a = binom.cdf(k1_low, cov1, f_bar) * binom.cdf(k2_high_ref_interpretation, cov2, 1 - f_bar)

    k1_high_ref_interpretation = np.floor( (1.0 - freq_thresh_high) * cov1) # = floor(0.3 * cov1)
    k2_low = np.floor(freq_thresh_low * cov2)
    prob_b = binom.cdf(k1_high_ref_interpretation, cov1, 1 - f_bar) * binom.cdf(k2_low, cov2, f_bar)

    return prob_a + prob_b

# --- Main Execution ---
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Starting global sweep analysis for all colonies...")

    # MODIFIED: Glob to find all colony SNV files
    snv_files = glob.glob(os.path.join(SNV_DIR, "snv_frequency_trajectory_*.csv"))
    if not snv_files:
        print(f"Error: No SNV files found in {SNV_DIR}")
        exit()

    print(f"Found {len(snv_files)} colony files to process.")

    all_interval_results = [] # Global list to store results from all colonies

    # --- MODIFIED: Outer loop over all colony files ---
    for snv_file in tqdm(snv_files, desc="Processing Colonies"):
        colony_id = os.path.basename(snv_file).split('_')[-1].replace('.csv', '')
        
        try:
            # --- Load Data ---
            dtype_spec = {
                'pos_in_contig': 'Int64', 'coverage': 'float64',
                'A': 'float64', 'T': 'float64', 'C': 'float64', 'G': 'float64',
                'departure_from_polarized_consensus': 'float64'
            }
            snv_df = pd.read_csv(snv_file, dtype=dtype_spec, low_memory=False)

            # --- Data Cleaning (same as before) ---
            snv_df['Month'] = pd.Categorical(snv_df['Month'], categories=MONTH_ORDER, ordered=True)
            snv_df['pos_in_contig'] = pd.to_numeric(snv_df['pos_in_contig'], errors='coerce')
            snv_df['coverage'] = pd.to_numeric(snv_df['coverage'], errors='coerce')
            snv_df['departure_from_polarized_consensus'] = pd.to_numeric(snv_df['departure_from_polarized_consensus'], errors='coerce')
            snv_df['alt_allele_count'] = (snv_df['departure_from_polarized_consensus'] * snv_df['coverage']).round().astype('Int64')
            snv_df['unique_SNV_identifier'] = snv_df['MAG_id'].astype(str) + '_' + \
                                              snv_df['contig_name'].astype(str) + '_' + \
                                              snv_df['pos_in_contig'].astype(str)
            essential_cols = ['MAG_id', 'Month', 'unique_SNV_identifier', 'coverage', 'departure_from_polarized_consensus', 'alt_allele_count']
            snv_df.dropna(subset=essential_cols, inplace=True)
            snv_df['coverage'] = snv_df['coverage'].astype(int)
            
            if snv_df.empty:
                tqdm.write(f"Colony {colony_id}: SKIPPED (empty after loading/cleaning)")
                continue

            # --- Pre-calculate Median Coverages (per colony) ---
            snv_df['median_mag_month_cov'] = snv_df.groupby(['MAG_id', 'Month'])['coverage'].transform('median')
            snv_df.dropna(subset=['median_mag_month_cov'], inplace=True)

            # --- Pre-filter 1: Sites always at low frequency (per colony) ---
            max_freq_per_site = snv_df.groupby('unique_SNV_identifier')['departure_from_polarized_consensus'].transform('max')
            snv_df_filtered = snv_df[max_freq_per_site > MIN_FREQ_THRESHOLD].copy()
            del snv_df # Free memory

            if snv_df_filtered.empty:
                tqdm.write(f"Colony {colony_id}: SKIPPED (empty after low-freq filter)")
                continue

            # --- Iterate through MAGs and Timepoint Pairs (per colony) ---
            all_mag_ids = snv_df_filtered['MAG_id'].unique()
            
            for mag_id in all_mag_ids: # No tqdm here, outer loop has it
                mag_df = snv_df_filtered[snv_df_filtered['MAG_id'] == mag_id].copy()
                mag_df.sort_values('Month', inplace=True)
                timepoints = mag_df['Month'].unique().tolist()

                for i in range(len(timepoints) - 1):
                    t1 = timepoints[i]
                    t2 = timepoints[i+1]

                    df_t1 = mag_df[mag_df['Month'] == t1].set_index('unique_SNV_identifier')
                    df_t2 = mag_df[mag_df['Month'] == t2].set_index('unique_SNV_identifier')
                    common_sites = df_t1.index.intersection(df_t2.index)
                    if len(common_sites) == 0: continue

                    df_t1_common = df_t1.loc[common_sites]
                    df_t2_common = df_t2.loc[common_sites]

                    # --- Pairwise Filters (same as before) ---
                    rel_cov_t1 = df_t1_common['coverage'] / df_t1_common['median_mag_month_cov'].replace(0, 1e-9)
                    rel_cov_t2 = df_t2_common['coverage'] / df_t2_common['median_mag_month_cov'].replace(0, 1e-9)
                    ratio_of_ratios = rel_cov_t2.replace(0, 1e-9) / rel_cov_t1.replace(0, 1e-9)
                    coverage_filter_mask = (ratio_of_ratios <= COVERAGE_RATIO_FOLD_CHANGE) & (ratio_of_ratios >= 1.0 / COVERAGE_RATIO_FOLD_CHANGE)
                    
                    pi_calc_df = pd.DataFrame({'cov1': df_t1_common['coverage'], 'cov2': df_t2_common['coverage'], 'alt1': df_t1_common['alt_allele_count'], 'alt2': df_t2_common['alt_allele_count']})
                    pi_values = pi_calc_df.apply(lambda row: calculate_pi(row['cov1'], row['cov2'], row['alt1'], row['alt2'], LOW_FREQ_SWEEP, HIGH_FREQ_SWEEP), axis=1)
                    pi_filter_mask = pi_values <= PI_THRESHOLD
                    
                    final_site_mask = coverage_filter_mask & pi_filter_mask
                    filtered_sites = common_sites[final_site_mask]
                    if len(filtered_sites) == 0: continue

                    df_t1_filt = df_t1_common.loc[filtered_sites]
                    df_t2_filt = df_t2_common.loc[filtered_sites]
                    pi_values_filt = pi_values.loc[filtered_sites]

                    # --- Calculate N_err and N_obs ---
                    n_err = pi_values_filt.sum()
                    case_a_mask = (df_t1_filt['departure_from_polarized_consensus'] <= LOW_FREQ_SWEEP) & (df_t2_filt['departure_from_polarized_consensus'] >= HIGH_FREQ_SWEEP)
                    case_b_mask = (df_t1_filt['departure_from_polarized_consensus'] >= HIGH_FREQ_SWEEP) & (df_t2_filt['departure_from_polarized_consensus'] <= LOW_FREQ_SWEEP)
                    observed_sweep_mask = case_a_mask | case_b_mask
                    n_obs = observed_sweep_mask.sum()
                    fdr_est = n_err / n_obs if n_obs > 0 else np.inf
                    poisson_p = poisson.sf(max(0, n_obs - 1), n_err) if n_err >= 0 else 1.0

                    # --- Store Results ---
                    all_interval_results.append({
                        'colony_id': colony_id, # MODIFIED: Store colony_id
                        'MAG_id': mag_id,
                        't1': t1,
                        't2': t2,
                        'n_sites_compared': len(filtered_sites),
                        'n_obs': n_obs,
                        'n_err': n_err,
                        'fdr_est': fdr_est,
                        'poisson_p': poisson_p
                    })
        except Exception as e:
            tqdm.write(f"Colony {colony_id}: FAILED during processing. Error: {e}")
            tqdm.write(traceback.format_exc())
            continue # Move to the next colony

    # --- MODIFIED: Global Post-processing (outside the colony loop) ---
    if not all_interval_results:
        print("\nNo timepoint comparisons yielded data across all colonies. Exiting.")
        exit()

    print("\nAll colonies processed. Performing global BH correction...")
    results_df = pd.DataFrame(all_interval_results)

    # Apply Benjamini-Hochberg correction to ALL Poisson p-values
    reject, pvals_corrected, _, _ = multipletests(results_df['poisson_p'].fillna(1.0), alpha=BH_ALPHA, method='fdr_bh')
    results_df['poisson_p_bh'] = pvals_corrected

    # Filter for significant intervals based on global results
    significant_intervals = results_df[
        (results_df['fdr_est'] < FDR_THRESHOLD) &
        (results_df['poisson_p_bh'] <= BH_ALPHA)
    ].copy()

    print(f"\nFound {len(significant_intervals)} globally significant MAG/interval pairs meeting criteria.")
    print(significant_intervals[['colony_id', 'MAG_id', 't1', 't2', 'n_obs', 'n_err', 'fdr_est', 'poisson_p_bh']])

    # --- MODIFIED: Identify SNVs involved in significant sweeps (re-loading files) ---
    print("\nIdentifying specific SNVs involved in significant shifts...")
    significant_snvs = []
    
    # Group by colony to load each file only once
    for colony_id, group in tqdm(significant_intervals.groupby('colony_id'), desc="Finding SNVs"):
        snv_file = os.path.join(SNV_DIR, f"snv_frequency_trajectory_{colony_id}.csv")
        try:
            # Load the corresponding colony file again
            snv_df = pd.read_csv(snv_file, dtype=dtype_spec, low_memory=False)
            # Re-run the same cleaning and pre-filtering
            snv_df['Month'] = pd.Categorical(snv_df['Month'], categories=MONTH_ORDER, ordered=True)
            snv_df['pos_in_contig'] = pd.to_numeric(snv_df['pos_in_contig'], errors='coerce')
            snv_df['coverage'] = pd.to_numeric(snv_df['coverage'], errors='coerce')
            snv_df['departure_from_polarized_consensus'] = pd.to_numeric(snv_df['departure_from_polarized_consensus'], errors='coerce')
            snv_df['alt_allele_count'] = (snv_df['departure_from_polarized_consensus'] * snv_df['coverage']).round().astype('Int64')
            snv_df['unique_SNV_identifier'] = snv_df['MAG_id'].astype(str) + '_' + snv_df['contig_name'].astype(str) + '_' + snv_df['pos_in_contig'].astype(str)
            snv_df.dropna(subset=essential_cols, inplace=True)
            snv_df['coverage'] = snv_df['coverage'].astype(int)
            snv_df['median_mag_month_cov'] = snv_df.groupby(['MAG_id', 'Month'])['coverage'].transform('median')
            snv_df.dropna(subset=['median_mag_month_cov'], inplace=True)
            max_freq_per_site = snv_df.groupby('unique_SNV_identifier')['departure_from_polarized_consensus'].transform('max')
            snv_df_filtered = snv_df[max_freq_per_site > MIN_FREQ_THRESHOLD].copy()
            del snv_df
            
            # Now iterate through the significant intervals for a colony
            for _, row in group.iterrows():
                mag_id = row['MAG_id']; t1 = row['t1']; t2 = row['t2']
                
                # Re-apply filters for this specific interval
                mag_df = snv_df_filtered[snv_df_filtered['MAG_id'] == mag_id].copy()
                df_t1 = mag_df[mag_df['Month'] == t1].set_index('unique_SNV_identifier')
                df_t2 = mag_df[mag_df['Month'] == t2].set_index('unique_SNV_identifier')
                common_sites = df_t1.index.intersection(df_t2.index)
                if len(common_sites) == 0: continue
                df_t1_common = df_t1.loc[common_sites]; df_t2_common = df_t2.loc[common_sites]
                
                rel_cov_t1 = df_t1_common['coverage'] / df_t1_common['median_mag_month_cov'].replace(0, 1e-9)
                rel_cov_t2 = df_t2_common['coverage'] / df_t2_common['median_mag_month_cov'].replace(0, 1e-9)
                ratio_of_ratios = rel_cov_t2.replace(0, 1e-9) / rel_cov_t1.replace(0, 1e-9)
                coverage_filter_mask = (ratio_of_ratios <= COVERAGE_RATIO_FOLD_CHANGE) & (ratio_of_ratios >= 1.0 / COVERAGE_RATIO_FOLD_CHANGE)

                pi_calc_df = pd.DataFrame({'cov1': df_t1_common['coverage'], 'cov2': df_t2_common['coverage'], 'alt1': df_t1_common['alt_allele_count'], 'alt2': df_t2_common['alt_allele_count']})
                pi_values = pi_calc_df.apply(lambda r: calculate_pi(r['cov1'], r['cov2'], r['alt1'], r['alt2'], LOW_FREQ_SWEEP, HIGH_FREQ_SWEEP), axis=1)
                pi_filter_mask = pi_values <= PI_THRESHOLD
                final_site_mask = coverage_filter_mask & pi_filter_mask
                filtered_sites = common_sites[final_site_mask]
                if len(filtered_sites) == 0: continue
                
                df_t1_filt = df_t1_common.loc[filtered_sites]; df_t2_filt = df_t2_common.loc[filtered_sites]

                # Identify the sweeping SNVs
                case_a_mask = (df_t1_filt['departure_from_polarized_consensus'] <= LOW_FREQ_SWEEP) & (df_t2_filt['departure_from_polarized_consensus'] >= HIGH_FREQ_SWEEP)
                case_b_mask = (df_t1_filt['departure_from_polarized_consensus'] >= HIGH_FREQ_SWEEP) & (df_t2_filt['departure_from_polarized_consensus'] <= LOW_FREQ_SWEEP)
                observed_sweep_mask = case_a_mask | case_b_mask
                sweeping_snv_ids = filtered_sites[observed_sweep_mask]

                for snv_id in sweeping_snv_ids:
                     freq1 = df_t1_filt.loc[snv_id, 'departure_from_polarized_consensus']
                     freq2 = df_t2_filt.loc[snv_id, 'departure_from_polarized_consensus']
                     contig = df_t1_filt.loc[snv_id, 'contig_name']
                     pos = df_t1_filt.loc[snv_id, 'pos_in_contig']
                     significant_snvs.append({
                         'colony_id': colony_id, 'MAG_id': mag_id, 't1': t1, 't2': t2,
                         'unique_SNV_identifier': snv_id, 'contig': contig, 'position': pos,
                         'freq_t1': freq1, 'freq_t2': freq2,
                         'fdr_est_interval': row['fdr_est'], 'poisson_p_bh_interval': row['poisson_p_bh']
                     })
        except Exception as e_load:
            tqdm.write(f"SKIPPED finding SNVs for colony {colony_id}. Error: {e_load}")
            continue

    # --- Save Global Results ---
    results_df.to_csv(os.path.join(OUTPUT_DIR, "sweep_interval_stats_ALL.csv"), index=False)
    print(f"\nSaved global interval statistics to: {os.path.join(OUTPUT_DIR, 'sweep_interval_stats_ALL.csv')}")

    if significant_snvs:
        significant_snvs_df = pd.DataFrame(significant_snvs)
        significant_snvs_df.to_csv(os.path.join(OUTPUT_DIR, "significant_sweep_snvs_ALL.csv"), index=False)
        print(f"Saved details of {len(significant_snvs_df)} SNVs in significant intervals to: {os.path.join(OUTPUT_DIR, 'significant_sweep_snvs_ALL.csv')}")
    else:
        print("\nNo specific SNVs met the final significance criteria across all intervals.")

    print(f"\n✅ Global sweep analysis complete.")