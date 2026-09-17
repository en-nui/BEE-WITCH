#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# SCRIPT METADATA
# ==============================================================================
# AUTHOR: Robin Chauhan (with assistance from Gemini)
# DATE: 2025-10-16
#
# PURPOSE:
# To create a master "gene map" that reconciles gene calls from anvi'o and DRAM,
# calculates detailed overlap metrics, and quantifies the number of "orphan"
# genes that are unique to each annotation set.
#
# WORKFLOW:
# 1.  LOAD DATA: It loads all required annotation and metadata files.
#
# 2.  CALCULATE OVERLAPS: Using an efficient Interval Tree approach, it finds the
#     best-matching DRAM gene for each anvi'o gene based on coordinate overlap.
#
# 3.  CALCULATE METRICS: It calculates detailed overlap statistics for each
#     reconciled gene pair.
#
# 4.  QUANTIFY ORPHANS: After reconciliation, it compares the set of matched
#     genes against the total set of genes from each caller to count the number
#     and proportion of genes that were not successfully matched.
#
# 5.  SAVE MAP & REPORT: It saves the final, detailed reconciliation table and
#     prints a summary of the orphan gene analysis.
#
# --- INPUTS ---
#   - RENAMED_CONTIG_TAXONOMY...csv: Master file for contig-to-MAG mapping.
#   - DRAM annotation files: The 'annotations.tsv' files from each colony.
#   - Anvi'o Universal Gene Calls: The master gene table from anvi'o.
#
# --- OUTPUTS ---
#   - gene_map_reconciliation.parquet: A file containing the detailed mapping.
#   - A summary printed to the console detailing the number of orphan genes.
# ==============================================================================

import pandas as pd
import os
import glob
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

# --- Configuration ---
CONTIG_TAXONOMY_FILE = "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"
DRAM_ANNOTATIONS_BASE_DIR = "/home/robinch/projects/BEE-WITCH/19_mMAG_annotations/"
ANVIO_GENES_FILE = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/scripts/diversity_parser/all_anvio_gene_calls_universal.txt"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/"
OUTPUT_GENE_MAP_FILE = os.path.join(OUTPUT_DIR, "gene_map_reconciliation.parquet")

# --- Main Execution ---
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # --- 1. Load All Necessary Data ---
    print("Loading all required input files...")
    contig_df = pd.read_csv(CONTIG_TAXONOMY_FILE, usecols=['updated_contig_id', 'MAG_id'])
    contig_to_mag_map = contig_df.set_index('updated_contig_id')['MAG_id'].to_dict()

    dram_files = glob.glob(os.path.join(DRAM_ANNOTATIONS_BASE_DIR, "colony_*", "DRAM_output", "annotations.tsv"))
    all_dram_genes = []
    for f in dram_files:
        try:
            dram_df = pd.read_csv(f, sep='\t', usecols=['scaffold', 'start_position', 'end_position'])
            dram_df.rename(columns={'scaffold': 'contig_name'}, inplace=True)
            all_dram_genes.append(dram_df)
        except Exception:
            continue
    dram_genes_df = pd.concat(all_dram_genes, ignore_index=True).drop_duplicates()
    dram_genes_df['unique_gene_id_DRAM'] = dram_genes_df['contig_name'] + '_' + dram_genes_df['start_position'].astype(str) + '_' + dram_genes_df['end_position'].astype(str)
    total_dram_genes = dram_genes_df['unique_gene_id_DRAM'].nunique()
    print(f"-> Found {total_dram_genes:,} unique DRAM gene calls.")

    anvio_genes_df = pd.read_csv(ANVIO_GENES_FILE, sep='\t', usecols=['contig', 'start', 'stop', 'unique_gene_callers_id'])
    anvio_genes_df.rename(columns={'contig': 'contig_name', 'start': 'start_position_anvio', 'stop': 'end_position_anvio'}, inplace=True)
    total_anvio_genes = anvio_genes_df['unique_gene_callers_id'].nunique()
    print(f"-> Found {total_anvio_genes:,} anvi'o gene calls.")

    # --- 2. Calculate Overlaps and Reconciliation Metrics ---
    print("\nCalculating overlaps and reconciliation metrics...")
    reconciled_results = []
    
    all_genes_grouped = pd.concat([dram_genes_df, anvio_genes_df]).groupby('contig_name')

    for contig, group in tqdm(all_genes_grouped, desc="Reconciling contigs"):
        dram_genes_on_contig = group.dropna(subset=['unique_gene_id_DRAM'])
        anvio_genes_on_contig = group.dropna(subset=['unique_gene_callers_id'])

        if dram_genes_on_contig.empty or anvio_genes_on_contig.empty:
            continue

        dram_tree = IntervalTree()
        for _, row in dram_genes_on_contig.iterrows():
            dram_tree.add(Interval(row['start_position'], row['end_position'], row.to_dict()))

        for _, anvio_row in anvio_genes_on_contig.iterrows():
            anvio_start, anvio_end = anvio_row['start_position_anvio'], anvio_row['end_position_anvio']
            overlapping_intervals = dram_tree.overlap(anvio_start, anvio_end)
            
            best_match_interval = None
            max_overlap = -1

            for interval in overlapping_intervals:
                overlap_size = min(anvio_end, interval.end) - max(anvio_start, interval.begin)
                if overlap_size > max_overlap:
                    max_overlap = overlap_size
                    best_match_interval = interval
            
            if best_match_interval:
                dram_match_data = best_match_interval.data
                anvio_range = anvio_end - anvio_start
                dram_range = dram_match_data['end_position'] - dram_match_data['start_position']
                
                if anvio_range > 0 and dram_range > 0:
                    percent_overlap_anvio = (max_overlap / anvio_range) * 100
                    percent_overlap_dram = (max_overlap / dram_range) * 100
                    average_overlap = (percent_overlap_anvio + percent_overlap_dram) / 2
                else:
                    average_overlap = 0.0

                reconciled_results.append({
                    'unique_gene_callers_id': anvio_row['unique_gene_callers_id'],
                    'unique_gene_id_DRAM': dram_match_data['unique_gene_id_DRAM'],
                    'anvio_range': anvio_range,
                    'dram_range': dram_range,
                    'overlap_bp': max_overlap,
                    'average_overlap_percent': average_overlap
                })

    # --- 3. Create Final Map and Save ---
    reconciled_df = pd.DataFrame(reconciled_results)
    
    reconciled_df_merged = pd.merge(reconciled_df, anvio_genes_df[['unique_gene_callers_id', 'contig_name']], on='unique_gene_callers_id')
    reconciled_df_merged['MAG_id'] = reconciled_df_merged['contig_name'].map(contig_to_mag_map)
    reconciled_df_merged.drop(columns=['contig_name'], inplace=True)
    
    reconciled_df_merged.to_parquet(OUTPUT_GENE_MAP_FILE)
    print(f"\n✅ Gene map created successfully with {len(reconciled_df):,} reconciled gene pairs.")
    print(f"Saved to: {OUTPUT_GENE_MAP_FILE}")

    # --- 4. NEW: Calculate and Report Orphan Gene Statistics ---
    reconciled_anvio_genes = reconciled_df['unique_gene_callers_id'].nunique()
    reconciled_dram_genes = reconciled_df['unique_gene_id_DRAM'].nunique()

    anvio_orphans = total_anvio_genes - reconciled_anvio_genes
    dram_orphans = total_dram_genes - reconciled_dram_genes

    anvio_orphan_prop = (anvio_orphans / total_anvio_genes) * 100 if total_anvio_genes > 0 else 0
    dram_orphan_prop = (dram_orphans / total_dram_genes) * 100 if total_dram_genes > 0 else 0

    print("\n--- Orphan Gene Analysis ---")
    print(f"Anvi'o Genes Not Found in DRAM: {anvio_orphans:,} ({anvio_orphan_prop:.2f}%)")
    print(f"DRAM Genes Not Found in Anvi'o: {dram_orphans:,} ({dram_orphan_prop:.2f}%)")