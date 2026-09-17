#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
from tqdm import tqdm

# ==============================================================================
# --- 1. CONFIGURATION ---
# ==============================================================================
CNV_RESULTS_PATH = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/anvio_based_cnv_outputs_numpy/all_mags_gene_coverage_stats_ANVIO_numpy.csv.gz"
METADATA_FILE_PATH = "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/all_filtered_mmag_contigs_concatenated.csv"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/pangenome_final_filtered/"
OUTPUT_CORE_GENES_PATH = os.path.join(OUTPUT_DIR, 'mag_colony_core_genes.csv.gz')
OUTPUT_ACCESSORY_GENES_PATH = os.path.join(OUTPUT_DIR, 'mag_colony_accessory_genes.csv.gz')
OUTPUT_SUMMARY_PATH = os.path.join(OUTPUT_DIR, 'pangenome_summary.csv')
OUTPUT_STABLE_GENES_PATH = os.path.join(OUTPUT_DIR, 'anvio_genes_of_interest.txt')


# --- PARAMETERS ---
CHUNK_SIZE = 1_000_000
CNV_LOWER_BOUND = 0.3
CNV_UPPER_BOUND = 3.0
PANGENOME_THRESHOLD = 0.90

# ==============================================================================
# --- 2. SCRIPT EXECUTION ---
# ==============================================================================
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("--- Starting New Pangenome Analysis with Advanced Filtering ---")

    # --- Pre-computation: Load small metadata ---
    metadata_df = pd.read_csv(METADATA_FILE_PATH, usecols=['MAG_id', 'Colony', 'sample_name'])
    metadata_df.rename(columns={'MAG_id': 'mag_id'}, inplace=True)
    unique_mags_initial = set(metadata_df['mag_id'].unique())
    mag_to_colony_map = metadata_df[['mag_id', 'Colony']].drop_duplicates()

    # --- Pass 1: Identify invalid MAG-sample pairs ---
    print("\n--- Pass 1: Identifying samples where MAG coverage is zero ---")
    invalid_pairs = set()
    chunk_iterator_1 = pd.read_csv(CNV_RESULTS_PATH, chunksize=CHUNK_SIZE, usecols=['mag_id', 'sample_id', 'mag_median_coverage'])
    for chunk in tqdm(chunk_iterator_1, desc="Scanning for invalid pairs"):
        invalid_chunk = chunk[chunk['mag_median_coverage'] == 0]
        if not invalid_chunk.empty:
            pairs = set(zip(invalid_chunk['mag_id'], invalid_chunk['sample_id']))
            invalid_pairs.update(pairs)
    affected_mags = {mag for mag, sample in invalid_pairs}
    proportion_removed = len(affected_mags) / len(unique_mags_initial)
    print(f"\n✅ Found {len(invalid_pairs):,} invalid MAG-sample combinations to exclude.")
    print(f"   This affected {len(affected_mags)} unique MAGs ({proportion_removed:.2%}).")

    # --- Pass 2: Calculate pangenome statistics on valid data ---
    print("\n--- Pass 2: Calculating core gene statistics ---")
    all_gene_stats = []
    chunk_iterator_2 = pd.read_csv(CNV_RESULTS_PATH, chunksize=CHUNK_SIZE)
    for chunk in tqdm(chunk_iterator_2, desc="Calculating gene stats"):
        is_valid_pair = ~pd.MultiIndex.from_frame(chunk[['mag_id', 'sample_id']]).isin(invalid_pairs)
        valid_chunk = chunk[is_valid_pair].copy()
        if valid_chunk.empty: continue
        valid_chunk['is_good'] = valid_chunk['gene_copy_number_median'].between(CNV_LOWER_BOUND, CNV_UPPER_BOUND)
        chunk_with_colony = pd.merge(valid_chunk, mag_to_colony_map, on='mag_id', how='left')
        chunk_with_colony.dropna(subset=['Colony'], inplace=True)
        chunk_stats = chunk_with_colony.groupby(['mag_id', 'Colony', 'unique_gene_callers_id']).agg(
            good_count=('is_good', 'sum'),
            total_count=('is_good', 'count')
        )
        all_gene_stats.append(chunk_stats)

    final_stats = pd.concat(all_gene_stats).groupby(level=[0, 1, 2]).sum()
    final_stats['prevalence_of_good'] = final_stats['good_count'] / final_stats['total_count']
    final_stats['is_core'] = final_stats['prevalence_of_good'] >= PANGENOME_THRESHOLD
    core_genes_set = set(final_stats[final_stats['is_core']].index.get_level_values('unique_gene_callers_id'))
    print(f"✅ Identified {len(core_genes_set)} core genes across all MAG/Colony pangenomes.")
    # Also save the list of stable genes (all genes considered in this step)
    stable_genes_set = set(final_stats.index.get_level_values('unique_gene_callers_id'))
    pd.Series(list(stable_genes_set)).to_csv(OUTPUT_STABLE_GENES_PATH, index=False, header=False)


    # --- Pass 3: Write final core and accessory files ---
    print("\n--- Pass 3: Writing final core and accessory files ---")
    for path in [OUTPUT_CORE_GENES_PATH, OUTPUT_ACCESSORY_GENES_PATH]:
        if os.path.exists(path):
            os.remove(path)
    header = True
    chunk_iterator_3 = pd.read_csv(CNV_RESULTS_PATH, chunksize=CHUNK_SIZE)
    for chunk in tqdm(chunk_iterator_3, desc="Writing output files"):
        is_valid_pair = ~pd.MultiIndex.from_frame(chunk[['mag_id', 'sample_id']]).isin(invalid_pairs)
        valid_chunk = chunk[is_valid_pair]
        if valid_chunk.empty: continue
        core_chunk = valid_chunk[valid_chunk['unique_gene_callers_id'].isin(core_genes_set)]
        accessory_chunk = valid_chunk[~valid_chunk['unique_gene_callers_id'].isin(core_genes_set) & valid_chunk['unique_gene_callers_id'].isin(stable_genes_set)] # Ensure accessory genes are from the stable set
        if not core_chunk.empty:
            core_chunk.to_csv(OUTPUT_CORE_GENES_PATH, mode='a', header=header, index=False, compression='gzip')
        if not accessory_chunk.empty:
            accessory_chunk.to_csv(OUTPUT_ACCESSORY_GENES_PATH, mode='a', header=header, index=False, compression='gzip')
        header = False

    # --- Final Summary (with updated column names) ---
    summary = final_stats.groupby(level=['mag_id', 'Colony']).agg(
        # 'is_core' is a boolean, count() gives total, sum() gives True counts
        genes_evaluated=('is_core', 'count'),
        core_genes_kept=('is_core', 'sum')
    ).reset_index()
    
    # Calculate removed genes and rename columns
    summary['accessory_genes_removed'] = summary['genes_evaluated'] - summary['core_genes_kept']
    summary.rename(columns={'mag_id': 'MAG_id'}, inplace=True)
    
    # Reorder columns to match request
    summary = summary[['MAG_id', 'Colony', 'genes_evaluated', 'core_genes_kept', 'accessory_genes_removed']]
    
    summary.to_csv(OUTPUT_SUMMARY_PATH, index=False)

    print("\n--- Pangenome Analysis Complete ---")
    print(f"✅ Core gene data saved to: {OUTPUT_CORE_GENES_PATH}")
    print(f"✅ Accessory gene data saved to: {OUTPUT_ACCESSORY_GENES_PATH}")
    print(f"✅ Pangenome summary saved to: {OUTPUT_SUMMARY_PATH}")
