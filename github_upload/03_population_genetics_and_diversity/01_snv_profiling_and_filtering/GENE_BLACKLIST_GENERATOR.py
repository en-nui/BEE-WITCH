#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# SCRIPT METADATA
# (Header is unchanged)
# ==============================================================================

import pandas as pd
import os
from tqdm import tqdm

# --- Configuration ---
TAXONOMY_FILE = "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"
GENES_OF_INTEREST_FILE = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/pangenome_final_filtered/anvio_genes_of_interest.txt"
ANVIO_GENES_FILE = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/scripts/diversity_parser/all_anvio_gene_calls_universal.txt"
OUTPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/blacklist_analysis/"
OUTPUT_FASTA = os.path.join(OUTPUT_DIR, "dereplicated_genes_for_clustering.fasta")

# --- Main Execution ---
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- 1. Load and Prepare Metadata ---
    print("Loading taxonomy and gene lists...")
    
    # --- FIX: Add 'updated_contig_id' to the usecols list ---
    tax_df = pd.read_csv(TAXONOMY_FILE, usecols=['MAG_id', 'Genus', 'type', 'updated_contig_id'])
    
    mag_to_genus = tax_df[tax_df['type'] == 'mmag'].set_index('MAG_id')['Genus'].to_dict()

    genes_of_interest = pd.read_csv(GENES_OF_INTEREST_FILE, header=None, names=['corresponding_gene_call'])
    genes_of_interest_set = set(genes_of_interest['corresponding_gene_call'])

    # --- 2. Load and Filter Anvi'o Gene Calls ---
    print("Loading and filtering anvi'o gene sequences...")
    anvio_genes_df = pd.read_csv(ANVIO_GENES_FILE, sep='\t', usecols=['gene_callers_id', 'unique_gene_callers_id', 'aa_sequence', 'contig'])
    
    # This line will now work correctly
    contig_to_mag = tax_df.set_index('updated_contig_id')['MAG_id'].to_dict()
    anvio_genes_df['MAG_id'] = anvio_genes_df['contig'].map(contig_to_mag)
    
    # Filter for genes of interest and valid MAGs
    anvio_genes_df = anvio_genes_df[anvio_genes_df['gene_callers_id'].isin(genes_of_interest_set)]
    anvio_genes_df.dropna(subset=['MAG_id', 'aa_sequence'], inplace=True)

    # --- 3. Dereplicate by unique_gene_callers_id ---
    print("Dereplicating gene sequences...")
    dereplicated_df = anvio_genes_df.drop_duplicates(subset=['unique_gene_callers_id'])
    dereplicated_df['Genus'] = dereplicated_df['MAG_id'].map(mag_to_genus)
    dereplicated_df.dropna(subset=['Genus'], inplace=True)
    
    print(f"-> Prepared {len(dereplicated_df)} unique gene sequences for clustering.")

    # --- 4. Write to FASTA file ---
    print(f"Writing sequences to FASTA file: {OUTPUT_FASTA}")
    with open(OUTPUT_FASTA, 'w') as f_out:
        for _, row in tqdm(dereplicated_df.iterrows(), total=len(dereplicated_df), desc="Writing FASTA"):
            # Header format: >unique_gene_callers_id|Genus
            header = f">{row['unique_gene_callers_id']}|{row['Genus']}\n"
            sequence = f"{row['aa_sequence']}\n"
            f_out.write(header)
            f_out.write(sequence)
            
    print("\n✅ Script 1 complete. You are now ready for the clustering step.")
