#!/usr/bin/env python3
# SCRIPT: 1_prepare_anvio_genes.py

import pandas as pd
import os

ANVIO_GENES_FILE = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/scripts/diversity_parser/all_anvio_gene_calls_universal.txt"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "anvio_genes_master_list.parquet")

print("--- Preparing Master Anvi'o Gene List ---")
anvio_df = pd.read_csv(
    ANVIO_GENES_FILE, 
    sep='\t', 
    usecols=['contig', 'start', 'stop', 'gene_callers_id', 'unique_gene_callers_id']
)
anvio_df.rename(columns={
    'contig': 'contig_name', 
    'start': 'anvio_start', 
    'stop': 'anvio_stop',
    'gene_callers_id': 'corresponding_gene_call'
}, inplace=True)

anvio_df.to_parquet(OUTPUT_FILE)
print(f"✅ Anvi'o master list saved to: {OUTPUT_FILE}")
