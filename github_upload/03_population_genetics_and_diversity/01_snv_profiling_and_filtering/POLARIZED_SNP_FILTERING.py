#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# SCRIPT METADATA
# ==============================================================================
# AUTHOR: Robin Chauhan (with assistance from Gemini)
# DATE: 2025-10-16
#
# SCRIPT 1 of 3: Initial SNP Cleaning (Memory-Optimized)
#
# PURPOSE:
# This script performs the initial cleaning of polarized SNP data. This version
# is optimized to handle extremely large files with minimal memory usage by
# avoiding the creation of large intermediate DataFrames.
#
# WORKFLOW (MEMORY-EFFICIENT MULTI-PASS):
# 1.  PASS 1 (CALCULATE 'D'): It scans all input files in chunks. Instead of
#     concatenating a large DataFrame, it collects raw coverage values for each
#     MAG/Month/Colony group into a dictionary. It then calculates the median
#     ('D') for each group directly from these lists of numbers.
#
# 2.  PASS 2 (FILTER & WRITE): It reads each input file a second time,
#     again in chunks. For each chunk, it filters SNPs and adds the
#     'unique_gene_callers_id', then appends the clean data directly to the
#     output file.
# ==============================================================================

import pandas as pd
import os
import glob
from tqdm import tqdm
from collections import defaultdict
import numpy as np

# --- CONFIGURATION ---
POLARIZED_SNPS_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/polarized_snps_final/"
ANVIO_GENES_PATH = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/scripts/diversity_parser/all_anvio_gene_calls_universal.txt"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/1_cleaned_snps/"
OUTPUT_SUMMARY_PATH = os.path.join(OUTPUT_DIR, 'cleaning_summary.csv')
CHUNK_SIZE = 500_000
SNP_COVERAGE_LOWER_FACTOR = 0.3
SNP_COVERAGE_UPPER_FACTOR = 3.0

# --- SCRIPT EXECUTION ---
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("--- Starting Script 1: Initial SNP Cleaning (Memory-Optimized) ---")
    all_files = sorted(glob.glob(os.path.join(POLARIZED_SNPS_DIR, '*.txt')))

    # --- Pre-computation: Create a fast lookup map for unique_gene_callers_id ---
    print("\nPre-loading anvi'o gene data to create a gene ID map...")
    try:
        anvio_genes_df = pd.read_csv(ANVIO_GENES_PATH, sep='\t', usecols=['gene_callers_id', 'unique_gene_callers_id'])
        gene_id_map = anvio_genes_df.set_index('gene_callers_id')['unique_gene_callers_id'].to_dict()
        print(f"✅ Created lookup map for {len(gene_id_map)} genes.")
    except FileNotFoundError:
        print(f"ERROR: Anvi'o gene file not found at '{ANVIO_GENES_PATH}'. Exiting."); exit()

    # --- Pass 1: Calculate 'D' for each MAG/Month/Colony (Memory-Safe Method) ---
    print("\n--- Pass 1: Calculating median coverage ('D') values ---")
    
    # Use a dictionary to collect coverage values, which is far more memory-efficient
    group_coverages = defaultdict(list)
    
    for f in tqdm(all_files, desc="Scanning for 'D'"):
        colony_id = os.path.basename(f).split('_')[1]
        chunk_iterator = pd.read_csv(f, sep='\t', chunksize=CHUNK_SIZE, low_memory=False, usecols=['MAG_id', 'sample_id', 'coverage', 'corresponding_gene_call'])
        for chunk in chunk_iterator:
            genic_chunk = chunk[chunk['corresponding_gene_call'] != -1]
            if genic_chunk.empty: continue
            
            genic_chunk['Colony'] = colony_id
            genic_chunk['Month'] = genic_chunk['sample_id'].str.split('_').str[1]
            
            # Append coverage values to lists instead of building a big DataFrame
            for mag, month, colony, cov in zip(genic_chunk['MAG_id'], genic_chunk['Month'], genic_chunk['Colony'], genic_chunk['coverage']):
                group_coverages[(mag, month, colony)].append(cov)

    if not group_coverages:
        print("ERROR: No genic SNPs found. Exiting."); exit()

    # Now, calculate the median for each group
    d_values_list = []
    for (mag, month, colony), cov_list in tqdm(group_coverages.items(), desc="Calculating Medians"):
        d_values_list.append({
            'MAG_id': mag,
            'Month': month,
            'Colony': colony,
            'D': np.median(cov_list)
        })
        
    d_values = pd.DataFrame(d_values_list)
    print(f"✅ Calculated 'D' values for {len(d_values)} groups.")

    # --- Pass 2: Filter and Write Clean SNP Files ---
    # This part of the logic is already memory-efficient and remains the same.
    print("\n--- Pass 2: Filtering SNPs and writing clean files ---")
    summary_stats = []
    for f in tqdm(all_files, desc="Cleaning files"):
        colony_id = os.path.basename(f).split('_')[1]
        filename = os.path.basename(f)
        output_file = os.path.join(OUTPUT_DIR, filename)
        if os.path.exists(output_file): os.remove(output_file)
        
        total_snps = 0
        removed_intergenic = 0
        removed_bad_coverage = 0
        final_kept = 0
        header = True
        
        chunk_iterator = pd.read_csv(f, sep='\t', chunksize=CHUNK_SIZE, low_memory=False)
        for chunk in chunk_iterator:
            total_snps += len(chunk)
            
            genic_mask = chunk['corresponding_gene_call'] != -1
            removed_intergenic += (~genic_mask).sum()
            genic_chunk = chunk[genic_mask].copy()
            if genic_chunk.empty: continue
            
            genic_chunk['unique_gene_callers_id'] = genic_chunk['corresponding_gene_call'].map(gene_id_map)
            genic_chunk['Month'] = genic_chunk['sample_id'].str.split('_').str[1]
            genic_chunk['Colony'] = colony_id
            chunk_with_d = pd.merge(genic_chunk, d_values, on=['MAG_id', 'Month', 'Colony'], how='left')
            
            min_cov = SNP_COVERAGE_LOWER_FACTOR * chunk_with_d['D']
            max_cov = SNP_COVERAGE_UPPER_FACTOR * chunk_with_d['D']
            coverage_mask = chunk_with_d['coverage'].between(min_cov, max_cov)
            removed_bad_coverage += (~coverage_mask).sum()
            
            clean_chunk = chunk_with_d[coverage_mask].drop(columns=['Month', 'Colony', 'D'])
            
            if not clean_chunk.empty:
                final_kept += len(clean_chunk)
                clean_chunk.to_csv(output_file, mode='a', header=header, sep='\t', index=False)
                header = False
                
        summary_stats.append({
            "colony_id": colony_id, "total_snps": total_snps, "removed_intergenic": removed_intergenic,
            "removed_bad_coverage": removed_bad_coverage, "final_snps_kept": final_kept
        })
        
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(OUTPUT_SUMMARY_PATH, index=False)
    
    print("\n✅ Script 1 complete.")
    print(f"Cleaned SNP files saved to: {OUTPUT_DIR}")
    print(f"Filtering summary saved to: {OUTPUT_SUMMARY_PATH}")