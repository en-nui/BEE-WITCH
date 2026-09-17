#!/usr/bin/env python3
# SCRIPT: 3_reconcile_gene_maps.py

import pandas as pd
import os
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

# --- Configuration ---
ANVIO_MASTER_LIST = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/anvio_genes_master_list.parquet"
DRAM_MASTER_LIST = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/dram_genes_master_list.parquet"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/"
OUTPUT_GENE_MAP_FILE = os.path.join(OUTPUT_DIR, "gene_map_reconciled.parquet")

print("--- Reconciling Anvi'o and DRAM Gene Maps ---")

anvio_df = pd.read_parquet(ANVIO_MASTER_LIST)
dram_df = pd.read_parquet(DRAM_MASTER_LIST)

reconciled_results = []

# Process one contig at a time to prevent any cross-contig errors
for contig, anvio_genes_on_contig in tqdm(anvio_df.groupby('contig_name'), desc="Reconciling contigs"):
    dram_genes_on_contig = dram_df[dram_df['contig_name'] == contig]
    
    if dram_genes_on_contig.empty:
        continue

    # Build an Interval Tree for the DRAM genes on this specific contig
    dram_tree = IntervalTree()
    for _, row in dram_genes_on_contig.iterrows():
        dram_tree.add(Interval(row['start_position'], row['end_position'], row.to_dict()))

    # Find the best match for each Anvi'o gene on this contig
    for _, anvio_row in anvio_genes_on_contig.iterrows():
        anvio_start, anvio_stop = anvio_row['anvio_start'], anvio_row['anvio_stop']
        overlapping_intervals = dram_tree.overlap(anvio_start, anvio_stop)
        
        best_match_interval = None
        max_overlap = -1

        for interval in overlapping_intervals:
            overlap_size = min(anvio_stop, interval.end) - max(anvio_start, interval.begin)
            if overlap_size > max_overlap:
                max_overlap = overlap_size
                best_match_interval = interval
        
        if best_match_interval:
            dram_match_data = best_match_interval.data
            
            # Create the final, merged record for this match
            result = {
                'corresponding_gene_call': anvio_row['corresponding_gene_call'],
                'unique_gene_callers_id': anvio_row['unique_gene_callers_id'],
                'contig_name': contig
            }
            # Add all columns from the matched DRAM gene
            result.update(dram_match_data)
            reconciled_results.append(result)

# Create the final dataframe and save
reconciled_df = pd.DataFrame(reconciled_results)
reconciled_df.to_parquet(OUTPUT_GENE_MAP_FILE)

print(f"\n✅ Gene map created successfully with {len(reconciled_df):,} reconciled gene pairs.")
print(f"Saved to: {OUTPUT_GENE_MAP_FILE}")
