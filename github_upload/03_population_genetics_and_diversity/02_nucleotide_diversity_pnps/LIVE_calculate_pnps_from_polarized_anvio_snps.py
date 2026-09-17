import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import product
import numpy as np
import os

# --- Helper Functions (Unchanged) ---
def parse_gene_calls(fasta_file_path, genes_to_keep=None):
    if not os.path.exists(fasta_file_path):
        print(f"Error: FASTA file not found at {fasta_file_path}")
        return None
    print(f"Reading and filtering gene sequences from: {fasta_file_path}")
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file_path, "fasta") if genes_to_keep is None or record.id in genes_to_keep}
    print(f"  -> Kept {len(sequences)} relevant gene sequences.")
    return sequences

def translate_codon(codon, is_first_codon=False):
    alt_starts = ["GTG", "TTG"]
    if is_first_codon and codon in alt_starts: return "M"
    try: return str(Seq(codon).translate())
    except: return "X"

def create_degeneracy_map():
    print("Pre-calculating degeneracy map...")
    degeneracy_map = {}
    bases = 'TCAG'
    codons = [''.join(p) for p in product(bases, repeat=3)]
    for codon in codons:
        try:
            original_aa = translate_codon(codon)
            if original_aa == '*': continue
        except: continue
        degeneracy_map[codon] = {}
        for pos in range(3):
            synonymous_changes = 0
            for new_base in bases:
                if new_base != codon[pos]:
                    new_codon_list = list(codon); new_codon_list[pos] = new_base
                    new_codon = "".join(new_codon_list)
                    if "*" in translate_codon(new_codon): continue
                    new_aa = translate_codon(new_codon)
                    if new_aa == original_aa: synonymous_changes += 1
            if synonymous_changes == 3: degeneracy_map[codon][pos + 1] = '4D'
            elif synonymous_changes == 1: degeneracy_map[codon][pos + 1] = '2D'
            else: degeneracy_map[codon][pos + 1] = '0D'
    return degeneracy_map

# --- Main Analysis Functions ---

def calculate_dn_ds_numerator(snp_df, gene_sequences_map):
    """Calculates dN and dS counts and returns a detailed log of all changes."""
    results = {gene_id: {'dN': 0, 'dS': 0} for gene_id in gene_sequences_map.keys()}
    detailed_changes = [] # NEW: List to store detailed per-SNP info

    print("Starting SNP processing loop for dN/dS counts...")
    for index, snp in snp_df.iterrows():
        gene_id_str = str(snp['corresponding_gene_call'])
        codon_order = snp['codon_order_in_gene']
        base_pos = snp['base_pos_in_codon']
        major_allele = snp['polarized_major_allele']
        minor_allele = snp['polarized_minor_allele']

        if pd.isna(major_allele) or pd.isna(minor_allele): continue
        
        gene_seq = gene_sequences_map[gene_id_str]
        codon_start_index = (codon_order - 1) * 3
        if codon_start_index < 0 or codon_start_index + 3 > len(gene_seq): continue
        
        template_codon = list(gene_seq[codon_start_index : codon_start_index + 3])
        
        ancestral_codon_list = template_codon[:]; ancestral_codon_list[base_pos - 1] = major_allele
        ancestral_codon = "".join(ancestral_codon_list)
        derived_codon_list = template_codon[:]; derived_codon_list[base_pos - 1] = minor_allele
        derived_codon = "".join(derived_codon_list)

        is_first = (codon_order == 1)
        ancestral_aa = translate_codon(ancestral_codon, is_first_codon=is_first)
        derived_aa = translate_codon(derived_codon, is_first_codon=False)
        
        if ancestral_aa == '*' or derived_aa == '*' or ancestral_aa == 'X' or derived_aa == 'X': continue

        change_type = 'S' if ancestral_aa == derived_aa else 'N'
        if change_type == 'S':
            results[gene_id_str]['dS'] += 1
        else:
            results[gene_id_str]['dN'] += 1
        
        # NEW: Append the detailed information to our log
        detailed_changes.append({
            'gene_id': gene_id_str,
            'codon_order': codon_order,
            'ancestral_codon': ancestral_codon,
            'derived_codon': derived_codon,
            'ancestral_aa': ancestral_aa,
            'derived_aa': derived_aa,
            'change_type': change_type
        })
    
    summary_df = pd.DataFrame.from_dict(results, orient='index').reset_index()
    summary_df.columns = ['corresponding_gene_call', 'dN', 'dS']
    summary_df['corresponding_gene_call'] = summary_df['corresponding_gene_call'].astype(int)
    
    detailed_df = pd.DataFrame(detailed_changes) # NEW: Create the detailed dataframe
    return summary_df, detailed_df

def calculate_final_ratios(dn_ds_df, denominator_file, genes_to_keep=None):
    # (This function is unchanged)
    print("\nLoading and filtering denominator data (N and S sites)...")
    try: sites_df = pd.read_csv(denominator_file, sep='\s+')
    except FileNotFoundError: print(f"Error: Denominator file not found at {denominator_file}"); return None
    if genes_to_keep:
        sites_df = sites_df[sites_df['corresponding_gene_call'].astype(str).isin(genes_to_keep)].copy()
        print(f"  -> Kept {len(sites_df)} relevant genes in denominator file.")
    final_df = pd.merge(dn_ds_df, sites_df, on='corresponding_gene_call', how='left')
    final_df.dropna(subset=['nN_gene_reference', 'nS_gene_reference'], inplace=True)
    final_df['pN'] = final_df['dN'] / final_df['nN_gene_reference']
    final_df['pS'] = final_df['dS'] / final_df['nS_gene_reference']
    with np.errstate(divide='ignore', invalid='ignore'):
        final_df['pN_pS_ratio'] = np.divide(final_df['pN'], final_df['pS'])
    final_df.replace([np.inf, -np.inf], 99, inplace=True)
    final_df.fillna(0, inplace=True)
    return final_df

# --- Main execution ---
if __name__ == '__main__':
    # 1. DEFINE YOUR FILE PATHS
    snp_table_file = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/annotated_with_unique_gene_id/colony_100_filtered_snps_cleaned_pooled.csv"
    gene_calls_fasta = "/home/robinch/projects/BEE-WITCH/40_internal_pnps_calculations_from_anvio/nt_gene_seq.fa"
    denominator_sites_file = "/home/robinch/projects/BEE-WITCH/22_anvio/Colony_100/100_merged_pnps_data/potentials.txt"
    output_snp_changes_file = "snp_codon_aa_changes.csv" # NEW: Output file for the detailed table

    # 2. PRE-FILTERING
    print(f"Loading SNP data to identify relevant genes from: {snp_table_file}")
    snp_df = pd.read_csv(snp_table_file)
    snp_df = snp_df[snp_df['in_coding_gene_call'].isin([1, True])].copy()
    genes_with_snps = set(snp_df['corresponding_gene_call'].astype(str).unique())
    print(f"Found {len(genes_with_snps)} genes with at least one coding SNP.")
    for col in ['codon_order_in_gene', 'base_pos_in_codon']:
        snp_df[col] = pd.to_numeric(snp_df[col], errors='coerce')
    snp_df.dropna(subset=['codon_order_in_gene', 'base_pos_in_codon'], inplace=True)
    snp_df = snp_df.astype({'codon_order_in_gene': int, 'base_pos_in_codon': int})

    # 3. RUN THE ANALYSIS
    gene_sequences = parse_gene_calls(gene_calls_fasta, genes_to_keep=genes_with_snps)
    if gene_sequences:
        dn_ds_counts, detailed_snp_report = calculate_dn_ds_numerator(snp_df, gene_sequences) # NEW: Unpack two dataframes
        if dn_ds_counts is not None:
            final_pn_ps_results = calculate_final_ratios(dn_ds_counts, denominator_sites_file, genes_to_keep=genes_with_snps)
            if final_pn_ps_results is not None:
                print("\n\n--- FINAL pN/pS RESULTS (First 20 Genes) ---")
                print(final_pn_ps_results.head(20).to_string())
                final_pn_ps_results.to_csv("final_pn_ps_ratios_colony_100.csv", index=False)
            
            # NEW: Save the detailed SNP report
            if detailed_snp_report is not None:
                print(f"\nSaving detailed SNP change report to: {output_snp_changes_file}")
                detailed_snp_report.to_csv(output_snp_changes_file, index=False)
