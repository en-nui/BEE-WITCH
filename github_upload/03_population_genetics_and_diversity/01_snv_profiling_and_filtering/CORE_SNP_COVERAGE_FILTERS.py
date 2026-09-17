#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# SCRIPT METADATA
# ==============================================================================
# AUTHOR: Robin Chauhan (with assistance from Gemini)
# DATE: 2025-10-15
#
# PURPOSE:
# This script performs a final, three-stage filtering on SNP datasets to produce
# a high-confidence set of variants for downstream population genetics. It is
# designed to be highly memory-efficient to handle massive input files.
#
# WORKFLOW (MULTI-PASS CHUNKING):
# 1.  PRE-COMPUTATION (Pass 1 over CNV data):
#     - It scans the large CNV results file to identify all MAG/sample pairs
#       that are unreliable (mag_median_coverage < 5). These are flagged for removal.
#     - It simultaneously builds a lookup map of the valid mag_median_coverage
#       for all other MAG/sample pairs.
#
# 2.  SINGLETON IDENTIFICATION (Pass 1 over SNP data):
#     - It reads each large colony SNP file in chunks.
#     - It applies the MAG quality and SNP coverage filters to each chunk.
#     - It counts the occurrences of all SNPs that pass these initial filters
#       to identify which ones are "singletons" (appearing in only one sample).
#
# 3.  FINAL FILTERING & WRITING (Pass 2 over SNP data):
#     - It reads each colony SNP file a final time, again in chunks.
#     - It applies all three filters concurrently:
#         a) Remove SNPs from low-quality MAG/sample pairs.
#         b) Remove SNPs with extreme coverage relative to the MAG median.
#         c) Remove singleton SNPs.
#     - Only SNPs that pass all three checks are written to the final output files.
#
# --- INPUTS ---
#   - CNV Results GZ: The large output file from the CNV calling script.
#   - Filtered SNP TXT files: The output from the previous script ('final_snps_by_colony/').
#
# --- OUTPUTS ---
#   - Final Polarized SNP TXT files: One per colony, containing the final, clean SNP set.
#   - Removal Summary CSV: A report detailing how many SNPs were removed at each stage.
# ==============================================================================

import pandas as pd
import os
from tqdm import tqdm
import glob
from collections import Counter

# ==============================================================================
# --- 1. CONFIGURATION ---
# ==============================================================================
# --- INPUT FILES ---
CNV_RESULTS_PATH = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/anvio_based_cnv_outputs_numpy/all_mags_gene_coverage_stats_ANVIO_numpy.csv.gz"
FILTERED_SNPS_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/final_snps_by_colony/"
METADATA_FILE_PATH = "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/all_filtered_mmag_contigs_concatenated.csv"

# --- OUTPUT FILES ---
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/polarized_snps_by_colony/"
OUTPUT_SUMMARY_PATH = os.path.join(OUTPUT_DIR, 'final_snp_filtering_summary.csv')

# --- PARAMETERS ---
CHUNK_SIZE = 500_000
MAG_COVERAGE_THRESHOLD = 5.0
SNP_COVERAGE_LOWER_FACTOR = 0.3
SNP_COVERAGE_UPPER_FACTOR = 3.0

# ==============================================================================
# --- 2. PRE-COMPUTATION: BUILD FILTERS (PASS 1 over CNV data) ---
# ==============================================================================
os.makedirs(OUTPUT_DIR, exist_ok=True)
print("--- Starting Final SNP Filtering Workflow ---")

print("\n--- Pass 1 (CNV Data): Building MAG/Sample filters ---")
invalid_pairs = set()
coverage_map = {}
chunk_iterator_1 = pd.read_csv(CNV_RESULTS_PATH, chunksize=CHUNK_SIZE, usecols=['mag_id', 'sample_id', 'mag_median_coverage'])

for chunk in tqdm(chunk_iterator_1, desc="Scanning CNV data"):
    # Identify pairs where MAG coverage is too low
    invalid_chunk = chunk[chunk['mag_median_coverage'] < MAG_COVERAGE_THRESHOLD]
    if not invalid_chunk.empty:
        pairs = set(zip(invalid_chunk['mag_id'], invalid_chunk['sample_id']))
        invalid_pairs.update(pairs)
    
    # Store the median coverage for all other valid pairs
    valid_chunk = chunk[chunk['mag_median_coverage'] >= MAG_COVERAGE_THRESHOLD].drop_duplicates(subset=['mag_id', 'sample_id'])
    for mag, sample, cov in zip(valid_chunk['mag_id'], valid_chunk['sample_id'], valid_chunk['mag_median_coverage']):
        coverage_map[(mag, sample)] = cov

print(f"✅ Flagged {len(invalid_pairs):,} MAG/sample pairs for removal due to low coverage.")

# --- Build contig-to-MAG map for linking SNPs to MAGs ---
print("\nBuilding contig-to-MAG map...")
metadata_df = pd.read_csv(METADATA_FILE_PATH, usecols=['MAG_id', 'updated_contig_id'])
contig_to_mag_map = metadata_df.set_index('updated_contig_id')['MAG_id'].to_dict()

# ==============================================================================
# --- 3. PROCESS SNP FILES (PASS 2 & 3 over SNP data) ---
# ==============================================================================
print("\n--- Processing SNP files for each colony ---")
snp_files = sorted(glob.glob(os.path.join(FILTERED_SNPS_DIR, '*_final_snps.txt')))
summary_stats = []

for snp_file in tqdm(snp_files, desc="Processing Colonies"):
    colony_id = os.path.basename(snp_file).replace('_final_snps.txt', '')
    output_file = os.path.join(OUTPUT_DIR, f'{colony_id}_polarized_snps.txt')

    total_snps = 0
    removed_by_mag_cov = 0
    removed_by_snp_cov = 0
    removed_by_singleton = 0
    
    # --- Pass 2 (SNP Data): Identify singletons ---
    snp_counter = Counter()
    chunk_iterator_2 = pd.read_csv(snp_file, sep='\t', chunksize=CHUNK_SIZE)
    for chunk in chunk_iterator_2:
        total_snps += len(chunk)
        chunk['mag_id'] = chunk['contig_name'].map(contig_to_mag_map)
        
        # Filter 1: MAG coverage
        chunk['pair'] = list(zip(chunk['mag_id'], chunk['sample_id']))
        mag_cov_mask = ~chunk['pair'].isin(invalid_pairs)
        chunk_after_mag_filter = chunk[mag_cov_mask].copy()
        
        # Filter 2: SNP coverage
        chunk_after_mag_filter['mag_median_coverage'] = chunk_after_mag_filter['pair'].map(coverage_map)
        min_cov = SNP_COVERAGE_LOWER_FACTOR * chunk_after_mag_filter['mag_median_coverage']
        max_cov = SNP_COVERAGE_UPPER_FACTOR * chunk_after_mag_filter['mag_median_coverage']
        snp_cov_mask = chunk_after_mag_filter['coverage'].between(min_cov, max_cov)
        
        passed_initial_filters = chunk_after_mag_filter[snp_cov_mask]
        
        # Count SNP occurrences
        snp_positions = zip(passed_initial_filters['contig_name'], passed_initial_filters['pos'])
        snp_counter.update(snp_positions)

    singletons = {snp for snp, count in snp_counter.items() if count == 1}

    # --- Pass 3 (SNP Data): Apply all filters and write ---
    if os.path.exists(output_file): os.remove(output_file)
    header = True
    final_kept_count = 0
    chunk_iterator_3 = pd.read_csv(snp_file, sep='\t', chunksize=CHUNK_SIZE)
    for chunk in chunk_iterator_3:
        chunk['mag_id'] = chunk['contig_name'].map(contig_to_mag_map)
        chunk['pair'] = list(zip(chunk['mag_id'], chunk['sample_id']))
        
        # Apply filters again to get masks
        mag_cov_mask = ~chunk['pair'].isin(invalid_pairs)
        
        chunk_after_mag_filter = chunk[mag_cov_mask].copy()
        chunk_after_mag_filter['mag_median_coverage'] = chunk_after_mag_filter['pair'].map(coverage_map)
        min_cov = SNP_COVERAGE_LOWER_FACTOR * chunk_after_mag_filter['mag_median_coverage']
        max_cov = SNP_COVERAGE_UPPER_FACTOR * chunk_after_mag_filter['mag_median_coverage']
        snp_cov_mask = chunk_after_mag_filter['coverage'].between(min_cov, max_cov)
        
        chunk_after_snp_filter = chunk_after_mag_filter[snp_cov_mask].copy()
        chunk_after_snp_filter['snp_pos'] = list(zip(chunk_after_snp_filter['contig_name'], chunk_after_snp_filter['pos']))
        singleton_mask = ~chunk_after_snp_filter['snp_pos'].isin(singletons)
        
        # Final clean chunk
        final_chunk = chunk_after_snp_filter[singleton_mask].drop(columns=['mag_id', 'pair', 'mag_median_coverage', 'snp_pos'])
        
        if not final_chunk.empty:
            final_chunk.to_csv(output_file, mode='a', header=header, sep='\t', index=False)
            header = False
        
        # Update counts for summary
        removed_by_mag_cov += (~mag_cov_mask).sum()
        removed_by_snp_cov += (mag_cov_mask & ~snp_cov_mask).sum()
        removed_by_singleton += (mag_cov_mask & snp_cov_mask & ~singleton_mask).sum()
        final_kept_count += len(final_chunk)
        
    summary_stats.append({
        "colony_id": colony_id,
        "total_snps_before_filtering": total_snps,
        "removed_low_mag_coverage": removed_by_mag_cov,
        "removed_bad_snp_coverage": removed_by_snp_cov,
        "removed_singleton_snps": removed_by_singleton,
        "final_snps_kept": final_kept_count
    })

# --- Save the final summary ---
summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv(OUTPUT_SUMMARY_PATH, index=False)

print("\n✅ Final filtering complete.")
print(f"Polarizable SNP files saved to: {OUTPUT_DIR}")
print(f"Summary of removed SNPs saved to: {OUTPUT_SUMMARY_PATH}")
