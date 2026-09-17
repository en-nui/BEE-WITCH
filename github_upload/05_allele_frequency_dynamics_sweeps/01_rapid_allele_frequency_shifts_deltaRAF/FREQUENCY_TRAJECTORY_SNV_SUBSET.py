#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# SCRIPT: subset_for_trajectory.py
# ==============================================================================

import pandas as pd
import os
import glob
from tqdm import tqdm
import concurrent.futures

# --- Configuration ---
FINAL_INPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/colony_polarized_outputs/pooled_by_month_global_consensus/"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories/"
TIMEPOINT_PRESENCE_THRESHOLD = 0.75 # Keep SNPs present in at least 75% of timepoints

# --- Worker Function ---
def create_trajectory_subset(colony_id):
    """
    Applies coverage and timepoint presence filters to create a subset of data
    for analyzing SNV frequency trajectories.
    """
    try:
        input_file = os.path.join(FINAL_INPUT_DIR, f"colony_{colony_id}_pooled_final_input.csv")
        output_file = os.path.join(OUTPUT_DIR, f"snv_frequency_trajectory_{colony_id}.csv")
        
        df = pd.read_csv(input_file)
        
        # --- Filter 1: Coverage Filter ---
        # For each MAG/Month, find the median SNP coverage.
        df['median_snp_coverage'] = df.groupby(['MAG_id', 'Month'])['coverage'].transform('median')
        
        # Filter out SNPs with coverage too low or too high relative to the median.
        df_filtered = df[
            (df['coverage'] >= df['median_snp_coverage'] / 3) &
            (df['coverage'] <= df['median_snp_coverage'] * 3)
        ].copy()
        
        # --- Filter 2: Timepoint Presence Filter ---
        # For each MAG, count the total number of unique months it appears in.
        mag_month_counts = df_filtered.groupby('MAG_id')['Month'].nunique().to_dict()
        df_filtered['total_months_for_mag'] = df_filtered['MAG_id'].map(mag_month_counts)
        
        # For each SNP site, count the number of months it's present.
        snp_site_month_counts = df_filtered.groupby(['MAG_id', 'contig_name', 'pos_in_contig'])['Month'].nunique().reset_index(name='months_present')
        
        # Merge this count back to the main dataframe.
        df_filtered = pd.merge(df_filtered, snp_site_month_counts, on=['MAG_id', 'contig_name', 'pos_in_contig'])
        
        # Keep only the SNP sites that meet the presence threshold.
        final_df = df_filtered[
            (df_filtered['months_present'] / df_filtered['total_months_for_mag']) >= TIMEPOINT_PRESENCE_THRESHOLD
        ].copy()

        # Clean up intermediate columns before saving.
        final_df.drop(columns=['median_snp_coverage', 'total_months_for_mag', 'months_present'], inplace=True)
        
        final_df.to_csv(output_file, index=False)
        return f"Colony {colony_id}: SUCCESS. Created trajectory subset."

    except Exception as e:
        return f"Colony {colony_id}: FAILED with error: {type(e).__name__} - {e}"

# --- Main Execution ---
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    all_files = glob.glob(os.path.join(FINAL_INPUT_DIR, "colony_*_pooled_final_input.csv"))
    colony_ids = sorted([os.path.basename(f).split('_')[1] for f in all_files])
    
    print(f"Found {len(colony_ids)} colonies to subset for trajectory analysis.")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        jobs = {executor.submit(create_trajectory_subset, cid): cid for cid in colony_ids}
        for future in tqdm(concurrent.futures.as_completed(jobs), total=len(colony_ids), desc="Creating Subsets"):
            print(f" -> {future.result()}")
            
    print("\n✅ Trajectory subsetting complete.")
