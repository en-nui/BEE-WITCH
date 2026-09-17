import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import glob
import concurrent.futures

# --- Configuration (Unchanged) ---
INPUT_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/annotated_with_unique_gene_id"
TAXONOMY_FILE = "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"
OUTPUT_DIR = os.path.join(os.path.dirname(INPUT_DIR), "popgen_summaries_final_filtered")
RECOVERY_THRESHOLD = 0.25

# --- Helper Functions (Unchanged) ---
def calculate_heterozygosity(df):
    df['heterozygosity'] = 0.0; valid_coverage = df['coverage'] > 1
    subset = df.loc[valid_coverage].copy()
    if not subset.empty:
        site_cov_adjustment = subset['coverage'] / (subset['coverage'] - 1)
        sum_of_squared_freqs = ((subset['A']/subset['coverage'])**2 + (subset['T']/subset['coverage'])**2 + (subset['C']/subset['coverage'])**2 + (subset['G']/subset['coverage'])**2)
        df.loc[valid_coverage, 'heterozygosity'] = (1 - sum_of_squared_freqs) * site_cov_adjustment
    return df

def calculate_popgen_stats(df, length_col, snv_col, coverage_col):
    df_out = df.copy(); n = df_out[coverage_col].astype(int).clip(lower=2); s = df_out[snv_col]
    max_n = n.max()
    if pd.isna(max_n) or max_n < 2: 
        df_out['Watterson_theta_per_bp'] = 0.0; df_out['TajimaD'] = np.nan
        return df_out
    i = np.arange(1, max_n); a1_series = np.cumsum(1 / i); a2_series = np.cumsum(1 / (i**2))
    a1 = a1_series[n - 2]; a2 = a2_series[n - 2]
    df_out['Watterson_theta_per_bp'] = np.divide(s / a1, df_out[length_col], where=(a1 > 0) & (df_out[length_col] > 0), out=np.zeros_like(s, dtype=float))
    pi = df_out['pi']; theta_w = s / a1
    b1 = (n + 1) / (3 * (n - 1)); b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1); c2 = b2 - (n + 2) / (a1 * n) + (a2 / a1**2)
    e1 = c1 / a1; e2 = c2 / (a1**2 + a2)
    var_d = (e1 * s) + (e2 * s * (s - 1))
    df_out['TajimaD'] = np.divide(pi - theta_w, np.sqrt(var_d), where=var_d > 0, out=np.full(len(df_out), np.nan))
    return df_out

# --- Worker Functions ---
def get_contig_info(file_path):
    try:
        df = pd.read_csv(file_path, usecols=['Colony', 'MAG_id', 'contig_name'])
        if df.empty: return None
        return df.drop_duplicates()
    except Exception: return None

def analyze_snp_data(file_path, whitelist_df):
    """Worker function for Pass 2, now with added safety check."""
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    df = pd.read_csv(file_path)
    df = pd.merge(df, whitelist_df, on=['Colony', 'MAG_id'], how='inner')
    if df.empty: return base_name, None, None

    if 'gene_length' in df.columns: df = df[df['gene_length'] >= 450].copy()
    if df.empty: return base_name, None, None
    
    grouping_keys_cov = ['MAG_id', 'Colony']
    df['median_coverage_group'] = df.groupby(grouping_keys_cov)['coverage'].transform('median')
    df = df[df['coverage'] <= (3 * df['median_coverage_group'])].copy()
    df = df[(df['departure_from_polarized_consensus'] >= 0.1) & (df['departure_from_polarized_consensus'] <= 0.9)]
    df.dropna(subset=['unique_gene_id'], inplace=True)
    if df.empty: return base_name, None, None
        
    df_with_het = calculate_heterozygosity(df)
    
    gene_level_groups = ['Colony', 'Month', 'MAG_id', 'unique_gene_id', 'gene_name', 'corresponding_gene_call']
    gene_df = df_with_het.groupby(gene_level_groups).agg(summed_heterozygosity=('heterozygosity','sum'),SNVs=('heterozygosity','size'),average_coverage=('coverage','mean'),gene_length=('gene_length','first')).reset_index()
    gene_df['pi'] = gene_df['summed_heterozygosity'] / gene_df['gene_length']
    gene_df['SNVs_per_100bp'] = (gene_df['SNVs'] / gene_df['gene_length']) * 100
    gene_df = calculate_popgen_stats(gene_df, 'gene_length', 'SNVs', 'average_coverage')
    
    # --- NEW SAFETY CHECK ---
    if gene_df.empty:
        # If no genes are left after aggregation, we can't create a MAG summary.
        # Return None for mag_df to prevent an error.
        return base_name, gene_df, None

    mag_level_groups = ['Colony', 'Month', 'MAG_id']
    mag_df = gene_df.groupby(mag_level_groups).agg(summed_heterozygosity=('summed_heterozygosity','sum'),total_SNVs=('SNVs','sum'),total_gene_length=('gene_length','sum'),average_coverage=('average_coverage','mean')).reset_index()
    mag_df['pi'] = mag_df['summed_heterozygosity'] / mag_df['total_gene_length']
    mag_df['SNVs_per_100bp'] = (mag_df['total_SNVs'] / mag_df['total_gene_length']) * 100
    mag_df = calculate_popgen_stats(mag_df, 'total_gene_length', 'total_SNVs', 'average_coverage')
    
    return base_name, gene_df, mag_df

# --- Main Execution Block (Unchanged) ---
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_files = glob.glob(os.path.join(INPUT_DIR, '*.csv'))
    print("--- Pass 1: Calculating Contig Recovery Rate ---")
    taxonomy_df = pd.read_csv(TAXONOMY_FILE)
    contigs_expected_df = taxonomy_df.groupby('MAG_id')['updated_contig_id'].nunique().reset_index()
    contigs_expected_df.rename(columns={'updated_contig_id': 'contigs_expected'}, inplace=True)
    all_contig_info = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        jobs = {executor.submit(get_contig_info, fp): fp for fp in all_files}
        for future in tqdm(concurrent.futures.as_completed(jobs), total=len(all_files), desc="Scanning contigs"):
            result = future.result()
            if result is not None: all_contig_info.append(result)
    if all_contig_info:
        observed_contigs_df = pd.concat(all_contig_info).drop_duplicates()
        contigs_observed_summary = observed_contigs_df.groupby(['Colony', 'MAG_id'])['contig_name'].nunique().reset_index()
        contigs_observed_summary.rename(columns={'contig_name': 'contigs_observed'}, inplace=True)
        recovery_summary_df = pd.merge(contigs_observed_summary, contigs_expected_df, on='MAG_id', how='left')
        recovery_summary_df['proportion_recovered'] = recovery_summary_df['contigs_observed'] / recovery_summary_df['contigs_expected']
        whitelist_df = recovery_summary_df[recovery_summary_df['proportion_recovered'] >= RECOVERY_THRESHOLD][['Colony', 'MAG_id']].copy()
        print(f"Whitelist created: {len(whitelist_df)} (Colony, MAG_id) pairs meet the {RECOVERY_THRESHOLD*100}% recovery threshold.")
        recovery_output_path = os.path.join(OUTPUT_DIR, "contig_recovery_summary.csv")
        recovery_summary_df.to_csv(recovery_output_path, index=False)
        print(f"Full contig recovery summary saved to: {recovery_output_path}")
    else:
        print("Could not generate contig recovery summary. Exiting.")
        exit()

    print("\n--- Pass 2: Running Population Genetics Analysis on Whitelisted Data ---")
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        jobs = {executor.submit(analyze_snp_data, fp, whitelist_df): fp for fp in all_files}
        for future in tqdm(concurrent.futures.as_completed(jobs), total=len(all_files), desc="Analyzing Files"):
            try:
                base_name, gene_summary_df, mag_summary_df = future.result()
                if gene_summary_df is not None:
                    gene_summary_df.to_csv(os.path.join(OUTPUT_DIR, f"{base_name}_by_gene.csv"), index=False)
                if mag_summary_df is not None:
                    mag_df.to_csv(os.path.join(OUTPUT_DIR, f"{base_name}_by_mag.csv"), index=False)
            except Exception as e:
                print(f"A file failed during analysis with error: {e}")
    print(f"\nProcessing complete. Final summary files are saved in: {OUTPUT_DIR}")