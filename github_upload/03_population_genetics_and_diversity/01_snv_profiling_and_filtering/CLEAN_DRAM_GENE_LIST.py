#!/usr/bin/env python3
# SCRIPT: 2_prepare_dram_genes.py

import pandas as pd
import os
import glob

DRAM_ANNOTATIONS_BASE_DIR = "/home/robinch/projects/BEE-WITCH/19_mMAG_annotations/"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "dram_genes_master_list.parquet")

print("--- Preparing Master DRAM Gene List with Annotations ---")
dram_cols_to_keep = [
    'scaffold', 'start_position', 'end_position', 
    'ko_id', 'pfam_hits', 'cazy_hits', 'vogdb_hits', 
    'viral_id', 'peptidase_hits', 'amr_hits', 'gene_id'
]

all_dram_genes = []
dram_files = glob.glob(os.path.join(DRAM_ANNOTATIONS_BASE_DIR, "colony_*", "DRAM_output", "annotations.tsv"))

for f in dram_files:
    try:
        temp_df = pd.read_csv(f, sep='\t')
        existing_cols = [col for col in dram_cols_to_keep if col in temp_df.columns]
        dram_df = temp_df[existing_cols]
        dram_df.rename(columns={'scaffold': 'contig_name'}, inplace=True)
        all_dram_genes.append(dram_df)
    except Exception:
        continue
        
dram_genes_df = pd.concat(all_dram_genes, ignore_index=True)
# Ensure we only have unique gene calls
dram_genes_df.drop_duplicates(subset=['contig_name', 'start_position', 'end_position'], inplace=True)

dram_genes_df.to_parquet(OUTPUT_FILE)
print(f"✅ DRAM master list saved to: {OUTPUT_FILE}")
