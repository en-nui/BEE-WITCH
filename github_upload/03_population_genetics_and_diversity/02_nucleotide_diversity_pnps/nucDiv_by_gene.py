import pandas as pd
import math
import os

################################################################################
# Module 1: Data Loading and Preparation
################################################################################

def load_and_prepare_data(metadata_path, variability_path):
    """
    Loads the metadata and variability data, filters metadata for vMAGs,
    and merges the two dataframes. Cleans up duplicated MAG_id columns.

    Args:
        metadata_path (str): Path to the metadata CSV file.
        variability_path (str): Path to the variability TXT file.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing the merged data with cleaned column names.
    """
    # Load the metadata
    metadata_df = pd.read_csv(metadata_path)

    # Filter metadata for viral MAGs (vMAGs)
    metadata_df_vMAG = metadata_df[metadata_df['type'] == "vMAG"].copy()

    # Load the Colony_100_NT variability data
    colony_df = pd.read_csv(variability_path, sep='\t')

    # Merge Colony_100_NT with metadata_df_vMAG
    merged_df = pd.merge(colony_df,
                           metadata_df_vMAG[['updated_contig_id', 'MAG_id', 'lifestyle', 'contig_length']], # Include 'lifestyle' and 'contig_length'
                           left_on='contig_name',
                           right_on='updated_contig_id',
                           how='left',
                           suffixes=('_x', '_y')) # Explicitly set suffixes

    # Remove MAG_id_y if it exists
    if 'MAG_id_y' in merged_df.columns:
        merged_df.drop(columns=['MAG_id_y'], inplace=True)

    # Rename MAG_id_x to MAG_id if it exists
    if 'MAG_id_x' in merged_df.columns:
        merged_df.rename(columns={'MAG_id_x': 'MAG_id'}, inplace=True)
    elif 'MAG_id' in merged_df.columns and 'MAG_id_x' not in merged_df.columns:
        # Handle case where merge might have resulted in just 'MAG_id' (e.g., if suffixes weren't needed)
        pass
    else:
        print("Warning: Could not find a MAG_id column after merging.")

    return merged_df

################################################################################
# Module 2: Heterozygosity Calculation
################################################################################

def calculate_heterozygosity(df):
    """
    Calculates per-site heterozygosity for a given DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame containing coverage and base counts (A, T, C, G).

    Returns:
        pandas.DataFrame: DataFrame with added 'site_cov', 'freq', and 'heterozygosity' columns.
    """
    # Calculate site coverage adjustment
    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)

    # Calculate the frequency component
    df['freq'] = df.apply(
        lambda row: (
            (row['A'] / row['coverage'])**2 +
            (row['T'] / row['coverage'])**2 +
            (row['C'] / row['coverage'])**2 +
            (row['G'] / row['coverage'])**2
        ),
        axis=1
    )

    # Calculate heterozygosity
    df['heterozygosity'] = (1 - df['freq']) * df['site_cov']

    return df

################################################################################
# Module 3: Execution and Filtering
################################################################################

if __name__ == "__main__":
    # Define the paths to your data files
    metadata_file_path = '/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_HOME_NO_TOUCH_20250327.csv'
    retained_file_path = '/home/robinch/projects/BEE-WITCH/21_final_r_vMAG_set/all_colonies_RETAINED_r_vMAGs.csv'
    anvio_base_dir = '/home/robinch/projects/BEE-WITCH/22_anvio/'
    output_nucdiv_dfs = {}
    filtering_summary_list = []

    # Iterate through subdirectories in the anvio base directory
    for subdir_name in os.listdir(anvio_base_dir):
        subdir_path = os.path.join(anvio_base_dir, subdir_name)
        if os.path.isdir(subdir_path):
            variability_file = None
            for filename in os.listdir(subdir_path):
                if filename.endswith("_annotated_variability.txt"):
                    variability_file = os.path.join(subdir_path, filename)
                    break  # Assuming only one variability file per subdir

            if variability_file:
                print(f"\n--- Processing variability file in subdirectory: {subdir_name} ---")
                print(f"  Variability file: {variability_file}")

                # Load and prepare the data
                merged_df_initial = load_and_prepare_data(metadata_file_path, variability_file)
                print(f"  Initial number of rows in merged_df: {len(merged_df_initial)}")
                initial_rows = len(merged_df_initial)

                # Filter out rows with NaN in gene_pos, start, or stop
                rows_before_gene_pos_nan_filter = len(merged_df_initial)
                merged_df_initial_no_gene_nan = merged_df_initial.dropna(subset=['gene_pos', 'start', 'stop']).copy()
                rows_after_gene_pos_nan_filter = len(merged_df_initial_no_gene_nan)
                gene_pos_nan_filtered = rows_before_gene_pos_nan_filter - rows_after_gene_pos_nan_filter
                print(f"  Number of rows after gene position NaN removal: {rows_after_gene_pos_nan_filter} (Filtered: {gene_pos_nan_filtered})")

                # Load the retained MAGs DataFrame
                retained_df = pd.read_csv(retained_file_path)
                merged_df_initial_no_gene_nan[['State', 'Month', 'Caste', 'Colony', 'Rep']] = merged_df_initial_no_gene_nan['sample_id'].str.split('_', expand=True)
                merged_df_initial_no_gene_nan['Colony'] = merged_df_initial_no_gene_nan['Colony'].astype(str)
                retained_df['Colony'] = retained_df['Colony'].astype(str)
                merged_df_initial_no_gene_nan = pd.merge(merged_df_initial_no_gene_nan, retained_df[['chrom', 'State', 'Month', 'Caste', 'Colony', 'Rep', 'mean']],
                                               left_on=['contig_name', 'State', 'Month', 'Caste', 'Colony', 'Rep'],
                                               right_on=['chrom', 'State', 'Month', 'Caste', 'Colony', 'Rep'],
                                               how='left')
                merged_df_initial_no_gene_nan.rename(columns={'mean': 'contig_mean'}, inplace=True)
                merged_df_initial_no_gene_nan.drop(columns=['chrom'], inplace=True)
                merged_df_initial_no_gene_nan['MAG_mean'] = merged_df_initial_no_gene_nan.groupby(['State', 'Month', 'Colony', 'Rep', 'MAG_id'])['contig_mean'].transform('mean')

                merged_df_filtered = merged_df_initial_no_gene_nan.copy() # Create a copy to perform filtering

                # --- Start of Simplified Filtering ---
                rows_before_contig_mean_nan_filter = len(merged_df_filtered)
                merged_df_filtered.dropna(subset=['contig_mean', 'MAG_mean'], inplace=True)
                rows_after_contig_mean_nan_filter = len(merged_df_filtered)
                contig_mean_nan_filtered = rows_before_contig_mean_nan_filter - rows_after_contig_mean_nan_filter
                print(f"  Number of rows after contig/MAG mean NaN removal: {rows_after_contig_mean_nan_filter} (Filtered: {contig_mean_nan_filtered})")

                filtered_df = merged_df_filtered # Assign the filtered DataFrame

                # --- End of Simplified Filtering ---

                # --- Start of Additional Coverage-Based Filtering ---
                initial_rows_coverage_filter = len(filtered_df)
                filtered_df = filtered_df[~((filtered_df['coverage'] > (filtered_df['contig_mean'] * 3)) | (filtered_df['coverage'] < (filtered_df['contig_mean'] * 0.3)))]
                rows_after_coverage_filter = len(filtered_df)
                coverage_filtered = initial_rows_coverage_filter - rows_after_coverage_filter
                print(f"  Number of rows after coverage-based filtering: {rows_after_coverage_filter} (Filtered: {coverage_filtered})")
                # --- End of Additional Coverage-Based Filtering ---

            # --- Start of Heterozygosity and Nucleotide Diversity Calculation at Gene Level ---
                df_with_heterozygosity_filtered = calculate_heterozygosity(filtered_df.copy()) # Calculate on a copy

                # --- Save the full filtered DataFrame ---
                output_dir = os.path.join(anvio_base_dir, subdir_name)
                full_filtered_filename = f"{subdir_name}_full_filtered.csv"
                full_filtered_filepath = os.path.join(output_dir, full_filtered_filename)
                try:
                    df_with_heterozygosity_filtered.to_csv(full_filtered_filepath, index=False)
                    print(f"  Full filtered DataFrame (per-site) written to: {full_filtered_filepath}")
                except Exception as e:
                    print(f"  Error writing full filtered DataFrame for subdirectory '{subdir_name}': {e}")

                nucdiv_by_gene_df = df_with_heterozygosity_filtered.groupby(['sample_id', 'MAG_id', 'updated_contig_id', 'gene_pos']).agg(
                    summed_heterozygosity=('heterozygosity', 'sum'),
                    SNVs=('heterozygosity', 'count'),
                    average_coverage=('coverage', 'mean'),
                    lifestyle=('lifestyle', 'first'), # Get the first lifestyle for each MAG
                    gene_start=('start', 'first'),
                    gene_stop=('stop', 'first')
                ).reset_index()

                # Calculate gene length
                nucdiv_by_gene_df['gene_length'] = nucdiv_by_gene_df['gene_stop'] - nucdiv_by_gene_df['gene_start']

                # Calculate nucleotide diversity
                nucdiv_by_gene_df['nucDiv_by_gene'] = nucdiv_by_gene_df['summed_heterozygosity'] / nucdiv_by_gene_df['gene_length']

                # Calculate number of SNPs per 100 bp
                nucdiv_by_gene_df['SNPs_per_100bp'] = nucdiv_by_gene_df['SNVs'] / (nucdiv_by_gene_df['gene_length'] / 100)

                # --- Calculate Watterson's Theta and Tajima's D ---
                def calculate_a1(n):
                    if n <= 1:
                        return 0
                    return sum(1/i for i in range(1, int(n)))
                def calculate_a2(n):
                    if n <= 1:
                        return 0
                    return sum(1/(i**2) for i in range(1, int(n)))
                def calculate_tajimas_d_var(n, s, a1, a2):
                    if n <= 1 or s <= 0:
                        return float('nan')
                    b1 = (n + 1) / (3 * (n - 1))
                    b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
                    c1 = b1 - 1 / a1
                    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1**2)
                    e1 = c1 / a1
                    e2 = c2 / (a1**2 + a2)
                    return (e1 * s) + (e2 * s * (s - 1))
                wtheta_list = []
                tajima_d_list = []
                for index, row in nucdiv_by_gene_df.iterrows():
                    n_snps = row['SNVs']
                    avg_coverage = row['average_coverage']
                    nuc_diversity = row['nucDiv_by_gene']
                    wtheta_denom = calculate_a1(avg_coverage)
                    wtheta = n_snps / wtheta_denom if wtheta_denom > 0 else 0
                    wtheta_list.append(wtheta)
                    a1 = calculate_a1(avg_coverage)
                    a2 = calculate_a2(avg_coverage)
                    tajima_d_var = calculate_tajimas_d_var(avg_coverage, n_snps, a1, a2)
                    tajima_d = (nuc_diversity - wtheta) / math.sqrt(tajima_d_var) if tajima_d_var > 0 else float('nan')
                    tajima_d_list.append(tajima_d)
                nucdiv_by_gene_df['Watterson_theta'] = wtheta_list
                nucdiv_by_gene_df['Wtheta_per_SNV'] = nucdiv_by_gene_df['Watterson_theta'] / nucdiv_by_gene_df['SNVs']
                nucdiv_by_gene_df['TajimaD'] = tajima_d_list

                print("  Nucleotide diversity by GENE with Watterson's Theta and Tajima's D:")
                print(nucdiv_by_gene_df[['sample_id', 'MAG_id', 'updated_contig_id', 'gene_pos', 'nucDiv_by_gene', 'SNPs_per_100bp', 'Watterson_theta', 'Wtheta_per_SNV', 'TajimaD']].head())

                # Store the resulting DataFrame
                output_nucdiv_dfs[subdir_name] = nucdiv_by_gene_df

                # Store filtering summary
                filtering_summary_list.append({
                    'Colony': subdir_name,
                    'Initial Rows': initial_rows,
                    'Gene Position NaN Filtered': gene_pos_nan_filtered,
                    'Contig/MAG Mean NaN Filtered': contig_mean_nan_filtered,
                    'Coverage Filtered': coverage_filtered,
                    'Final Rows (Gene Level)': len(nucdiv_by_gene_df),
                    'Rows After Full Filtering (Per-Site)': len(df_with_heterozygosity_filtered)
                })

            else:
                print(f"  Warning: No *_annotated_variability.txt file found in subdirectory: {subdir_name}")

    print("\n--- Summary of nucdiv_by_gene_df DataFrames created ---")
    for subdir, df in output_nucdiv_dfs.items():
        print(f"  Subdirectory: {subdir}, DataFrame shape: {df.shape}")

    # --- Write each DataFrame to its respective subdirectory ---
    print("\n--- Writing nucdiv_by_gene_df to respective subdirectories ---")
    for subdir_name, df in output_nucdiv_dfs.items():
        output_dir = os.path.join(anvio_base_dir, subdir_name)
        output_filename = f"{subdir_name}_nucdiv_by_gene.csv"
        output_filepath = os.path.join(output_dir, output_filename)
        try:
            df.to_csv(output_filepath, index=False)
            print(f"  Gene-level DataFrame for subdirectory '{subdir_name}' written to: {output_filepath}")
        except Exception as e:
            print(f"  Error writing gene-level DataFrame for subdirectory '{subdir_name}': {e}")

    print("\n--- Filtering Summary ---")
    filtering_summary_df = pd.DataFrame(filtering_summary_list)
    print(filtering_summary_df)

    print("\nFinished processing all subdirectories and writing output files.")
