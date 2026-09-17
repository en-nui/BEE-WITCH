#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# SCRIPT METADATA
# ==============================================================================
# AUTHOR: Robin Robinson
# DATE: 2025-10-15
#
# PURPOSE:
# This script performs a two-stage filtering process on large SNP datasets
# (from 'filtered_variability.txt' files) to retain only high-confidence SNPs
# that are suitable for population genetics analysis.
#
# WORKFLOW:
# 1.  IDENTIFY HIGH-QUALITY MAGs: It first reads a pangenome summary to identify
#     "high-quality" MAGs, defined as those having a sufficient number of core
#     genes (e.g., >= 100).
#
# 2.  BUILD CORE GENE REFERENCE: It then identifies the genomic coordinates of all
#     core genes that belong to these high-quality MAGs. This reference is
#     built into a highly efficient `IntervalTree` data structure for rapid
#     location lookups.
#
# 3.  FILTER SNPs IN CHUNKS: The script iterates through each colony's large SNP
#     file, reading it in memory-efficient chunks. For each chunk, it applies
#     two concurrent filters:
#       a) MAG Filter: It discards any SNP originating from a "low-quality" MAG.
#       b) Gene Filter: It discards any remaining SNP that does not fall within
#          the genomic coordinates of a core gene.
#
# 4.  GENERATE OUTPUTS: It writes the SNPs that pass both filters to new, clean
#     output files and generates a summary report detailing how many SNPs were
#     removed by each filter for every colony.
#
# KEY FEATURE - MEMORY EFFICIENCY:
# This script is specifically designed to handle massive input files (many GBs)
# without crashing. By processing the SNP data in small chunks, the RAM usage
# remains low and constant, making it suitable for large-scale analysis on
# standard computing hardware.
#
# --- INPUTS ---
#   - Pangenome Summary CSV: A summary file with MAG_id and core_genes_kept counts.
#   - Core Genes CSV: The detailed output from the pangenome script, listing all core genes.
#   - Anvi'o Gene Calls TXT: A master file with coordinates for all genes.
#   - Metadata CSV: A file mapping contig IDs to MAG IDs.
#   - Variability Files: The per-colony 'filtered_variability.txt' files containing raw SNP data.
#
# --- OUTPUTS ---
#   - Filtered SNP TXT files: A new set of SNP files, one per colony, containing only the high-confidence SNPs.
#   - Removal Summary CSV: A report quantifying the number of SNPs removed at each filtering stage per colony.
# ==============================================================================

import pandas as pd
import os
from intervaltree import Interval, IntervalTree
import glob
from tqdm import tqdm

# ==============================================================================
# --- 1. CONFIGURATION ---
# ==============================================================================
# --- INPUT FILES ---
PANGENOME_SUMMARY_PATH = '/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/pangenome_final_filtered/pangenome_summary.csv'
CORE_GENES_PATH = '/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/pangenome_final_filtered/mag_colony_core_genes.csv.gz'
ANVIO_GENES_PATH = '/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/scripts/diversity_parser/all_anvio_gene_calls_universal.txt'
METADATA_FILE_PATH = '/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/all_filtered_mmag_contigs_concatenated.csv'
VARIABILITY_BASE_DIR = '/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/'

# --- OUTPUT FILES ---
OUTPUT_DIR = '/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/final_snps_by_colony/'
OUTPUT_SUMMARY_PATH = os.path.join(OUTPUT_DIR, 'removed_snp_counts.csv')

# --- PARAMETERS ---
CORE_GENE_THRESHOLD = 100
CHUNK_SIZE = 500_000 # Process 500k SNPs at a time to keep memory low

# ==============================================================================
# --- 2. PRE-COMPUTATION: BUILD FILTERS ---
# ==============================================================================
os.makedirs(OUTPUT_DIR, exist_ok=True)
print("--- Starting SNP Filtering with Two Concurrent Filters ---")

# --- Filter 1: Identify "High-Quality" MAGs ---
print(f"\n--- Step 1: Identifying MAGs with >= {CORE_GENE_THRESHOLD} core genes ---")
try:
    pangenome_summary_df = pd.read_csv(PANGENOME_SUMMARY_PATH)
    high_core_mags = pangenome_summary_df[pangenome_summary_df['core_genes_kept'] >= CORE_GENE_THRESHOLD]
    high_core_mags_set = set(high_core_mags['MAG_id'].unique())
    print(f"Found {len(high_core_mags_set)} MAGs that pass the quality threshold.")
except FileNotFoundError as e:
    print(f"ERROR: Pangenome summary file not found. {e}"); exit()

# --- Filter 2: Build Interval Tree for Core Genes of High-Quality MAGs ---
print("\n--- Step 2: Building interval tree reference for core genes ---")
try:
    core_genes_df = pd.read_csv(CORE_GENES_PATH, usecols=['unique_gene_callers_id', 'mag_id'])
    core_genes_df = core_genes_df[core_genes_df['mag_id'].isin(high_core_mags_set)]
    core_gene_ids = set(core_genes_df['unique_gene_callers_id'].unique())
    
    all_genes_df = pd.read_csv(ANVIO_GENES_PATH, sep='\t')
    core_gene_coords = all_genes_df[all_genes_df['unique_gene_callers_id'].isin(core_gene_ids)]
except FileNotFoundError as e:
    print(f"ERROR: A required gene file was not found. {e}"); exit()

contig_interval_trees = {}
for contig, genes_on_contig in tqdm(core_gene_coords.groupby('contig'), desc="Building trees"):
    tree = IntervalTree()
    for _, gene in genes_on_contig.iterrows():
        tree.add(Interval(gene['start'], gene['stop']))
    contig_interval_trees[contig] = tree
print(f"Reference built for {len(contig_interval_trees)} contigs containing core genes.")

# --- Build a map from contig to MAG for efficient lookup ---
print("\nBuilding contig-to-MAG map...")
metadata_df = pd.read_csv(METADATA_FILE_PATH, usecols=['MAG_id', 'updated_contig_id'])
contig_to_mag_map = metadata_df.set_index('updated_contig_id')['MAG_id'].to_dict()

# ==============================================================================
# --- 3. PROCESS SNP FILES IN CHUNKS ---
# ==============================================================================
print("\n--- Step 3: Processing colony variability files in chunks ---")

def is_in_core_gene(row):
    """Helper function to query the interval tree for a single SNP."""
    contig = row['contig_name']
    pos = row['pos']
    if contig in contig_interval_trees:
        return bool(contig_interval_trees[contig].at(pos)) # Explicitly return boolean
    return False

colony_dirs = sorted(glob.glob(os.path.join(VARIABILITY_BASE_DIR, 'colony_*')))
summary_stats = []

for colony_dir in tqdm(colony_dirs, desc="Filtering Colonies"):
    colony_id = os.path.basename(colony_dir)
    variability_file = os.path.join(colony_dir, 'filtered_variability.txt')
    output_file = os.path.join(OUTPUT_DIR, f'{colony_id}_final_snps.txt')
    if not os.path.exists(variability_file): continue

    total_snps_in_file = 0
    removed_by_mag_filter = 0
    removed_by_gene_filter = 0
    final_snps_kept = 0
    header = True

    try:
        chunk_iterator = pd.read_csv(variability_file, sep='\t', chunksize=CHUNK_SIZE)
        for chunk in chunk_iterator:
            total_snps_in_file += len(chunk)
            
            chunk['mag_id'] = chunk['contig_name'].map(contig_to_mag_map)
            chunk.dropna(subset=['mag_id'], inplace=True)

            is_high_core_mag_mask = chunk['mag_id'].isin(high_core_mags_set)
            is_in_core_gene_mask = chunk.apply(is_in_core_gene, axis=1)
            
            removed_by_mag_filter += (~is_high_core_mag_mask).sum()
            removed_by_gene_filter += (is_high_core_mag_mask & ~is_in_core_gene_mask).sum()

            final_keep_mask = is_high_core_mag_mask & is_in_core_gene_mask
            filtered_chunk = chunk.loc[final_keep_mask, chunk.columns != 'mag_id'] # Drop temp mag_id col
            
            if not filtered_chunk.empty:
                final_snps_kept += len(filtered_chunk)
                filtered_chunk.to_csv(output_file, mode='a', header=header, sep='\t', index=False)
                header = False
    except Exception as e:
        print(f"ERROR processing {colony_id}: {e}. Skipping.")
        continue
    
    summary_stats.append({
        "colony_id": colony_id,
        "total_snps": total_snps_in_file,
        "snps_removed_low_core_mag": removed_by_mag_filter,
        "snps_removed_not_in_core_gene": removed_by_gene_filter,
        "final_snps_kept": final_snps_kept
    })

# --- Save the final summary of removed SNPs ---
summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv(OUTPUT_SUMMARY_PATH, index=False)

print("\n✅ All colonies processed successfully.")
print(f"Final SNP files saved to: {OUTPUT_DIR}")
print(f"Summary of removed SNPs saved to: {OUTPUT_SUMMARY_PATH}")
