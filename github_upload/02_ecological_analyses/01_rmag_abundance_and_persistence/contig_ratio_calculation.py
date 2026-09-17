import pandas as pd
import numpy as np
import os

# --- Configuration Paths ---
# !!! IMPORTANT: LOAD YOUR ACTUAL df_mMAG HERE !!!
# This script assumes df_mMAG is a pandas DataFrame available in your environment.
# Uncomment and modify the line below if you need to load it from a file:
# df_mMAG = pd.read_csv('/path/to/your/actual_df_mMAG.csv')

# --- IMPORTANT: Ensure df_mMAG is loaded or defined here ---
# If df_mMAG is not defined, this script will exit.
df_mMAG = pd.read_csv('/home/robinch/projects/BEE-WITCH/16_rMAG_mapping/colony_rMAG_mapping_with_metadata/concat_input/concat_annotated_summaries_for_R.csv')

try:
    # This check is just to provide a helpful message if df_mMAG isn't loaded.
    # In your actual use, ensure df_mMAG is loaded from your data source.
    if 'df_mMAG' not in locals() and 'df_mMAG' not in globals():
        raise NameError("df_mMAG DataFrame not found. Please load your main DataFrame (e.g., df_mMAG = pd.read_csv('your_file.csv')) before running this script.")
    
    print(f"Using df_mMAG with {df_mMAG.shape[0]} rows and {df_mMAG.shape[1]} columns.")
    
except NameError as e:
    print(f"Error: {e}")
    exit() # Exit if df_mMAG is not available

# Path to your RAW RENAMED_CONTIG_HOME file, crucial for metadata_sum_contig_length
RAW_RENAMED_FILE_PATH = '/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_HOME_NO_TOUCH_20250327.csv'

# Output directory for the new ratio table
OUTPUT_DIR = 'ratio_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("--- Starting Enhanced Ratio Analysis Script ---")

# Ensure 'contig_length' and 'mean' columns are numeric
df_mMAG['contig_length'] = pd.to_numeric(df_mMAG['contig_length'], errors='coerce')
df_mMAG['mean'] = pd.to_numeric(df_mMAG['mean'], errors='coerce')



## 1. Prepare Base DataFrame: `sample_id` Creation and `mMAG` Filtering

print("\n1. Creating 'sample_id' and filtering for 'mMAG' type...")
df_mMAG['sample_id'] = (df_mMAG['State'].astype(str).fillna('') + '_' +
                        df_mMAG['Month'].astype(str).fillna('') + '_' +
                        df_mMAG['Caste'].astype(str).fillna('') + '_' +
                        df_mMAG['Colony'].astype(str).fillna('') + '_' +
                        df_mMAG['Rep'].astype(str).fillna(''))
df_mMAG_filtered = df_mMAG[df_mMAG['type'] == 'mMAG'].copy()
print(f"   Base 'df_mMAG_filtered' created. Shape: {df_mMAG_filtered.shape}")

print("\n2. Aggregating summed_contig_length, average_mean, median_mean, and summed_mean per (sample_id, MAG_id, Genus)...")
df_aggregated_sample_mag_stats = df_mMAG_filtered.groupby(['sample_id', 'MAG_id', 'Genus']).agg(
    summed_contig_length_per_sample_mag=('contig_length', 'sum'),
    average_mean_per_sample_mag=('mean', 'mean'),  # Mean coverage
    median_mean_per_sample_mag=('mean', 'median'), # Median coverage
    summed_mean_per_sample_mag=('mean', 'sum')     # Summed mean coverage (new)
).reset_index()

# Drop rows where aggregated values might be NaN
df_aggregated_sample_mag_stats.dropna(
    subset=[
        'summed_contig_length_per_sample_mag',
        'average_mean_per_sample_mag',
        'median_mean_per_sample_mag',
        'summed_mean_per_sample_mag' # New column
    ],
    inplace=True
)
print(f"   'df_aggregated_sample_mag_stats' created. Shape: {df_aggregated_sample_mag_stats.shape}")
print("   Head of df_aggregated_sample_mag_stats:\n", df_aggregated_sample_mag_stats.head())

print(f"\n3. Dynamically calculating 'metadata_sum_contig_length' from raw file: {RAW_RENAMED_FILE_PATH}...")
df_metadata_contig_length_ref = pd.DataFrame() # Initialize as empty

try:
    df_raw_renamed = pd.read_csv(RAW_RENAMED_FILE_PATH)
    print(f"   Raw RENAMED file loaded. Shape: {df_raw_renamed.shape}")

    # Ensure 'contig_length' is numeric for summing
    df_raw_renamed['contig_length'] = pd.to_numeric(df_raw_renamed['contig_length'], errors='coerce')

    # Filter for mMAGs that are also representatives (MAG_id == r_MAG_id)
    df_mmag_representatives_ref = df_raw_renamed[
        (df_raw_renamed['type'] == 'mMAG') &
        (df_raw_renamed['MAG_id'] == df_raw_renamed['r_MAG_id'])
    ].copy()
    print(f"   Filtered raw data for mMAG representatives. Shape: {df_mmag_representatives_ref.shape}")

    # Group by r_MAG_id and sum their contig_lengths to get the metadata reference
    df_metadata_contig_length_ref = df_mmag_representatives_ref.groupby('r_MAG_id').agg(
        metadata_sum_contig_length=('contig_length', 'sum')
    ).reset_index()

    # Rename r_MAG_id to MAG_id for consistent merging key
    df_metadata_contig_length_ref.rename(columns={'r_MAG_id': 'MAG_id'}, inplace=True)

    # Drop rows where metadata_sum_contig_length might be NaN after aggregation
    df_metadata_contig_length_ref.dropna(subset=['metadata_sum_contig_length'], inplace=True)

    print(f"   'df_metadata_contig_length_ref' calculated. Shape: {df_metadata_contig_length_ref.shape}")
    print("   Head of df_metadata_contig_length_ref:\n", df_metadata_contig_length_ref.head())

except FileNotFoundError:
    print(f"   ERROR: Raw RENAMED file NOT FOUND at {RAW_RENAMED_FILE_PATH}. Cannot calculate metadata reference.")
    exit() # Exit if this critical file isn't found
except KeyError as e:
    print(f"   ERROR: Missing expected column '{e}' in raw RENAMED file. Cannot calculate metadata reference.")
    exit() # Exit if essential columns are missing
    
    print("\n4. Merging aggregated data and metadata reference to calculate ratios...")

# Merge df_aggregated_sample_mag_stats with df_metadata_contig_length_ref
# A left merge will keep all (sample_id, MAG_id, Genus) pairs and add the metadata sum if available.
df_ratio_analysis = df_aggregated_sample_mag_stats.merge(
    df_metadata_contig_length_ref,
    on='MAG_id',
    how='left'
)

# Fill NaN in 'metadata_sum_contig_length' with 0. This prepares for ratio calculations;
# ratios will be NaN if the denominator is 0.
df_ratio_analysis['metadata_sum_contig_length'].fillna(0, inplace=True)
print(f"   Merged table created. Shape: {df_ratio_analysis.shape}")
print("   Head of merged table (before ratio calculations):\n", df_ratio_analysis.head())

# Calculate contig_length_ratio, handling division by zero explicitly
df_ratio_analysis['contig_length_ratio'] = df_ratio_analysis.apply(
    lambda row: row['summed_contig_length_per_sample_mag'] / row['metadata_sum_contig_length']
                if row['metadata_sum_contig_length'] != 0 else np.nan, axis=1
)

# Calculate summed_mean_div_sample_contig_length, handling division by zero
df_ratio_analysis['summed_mean_div_sample_contig_length'] = df_ratio_analysis.apply(
    lambda row: row['summed_mean_per_sample_mag'] / row['summed_contig_length_per_sample_mag']
                if row['summed_contig_length_per_sample_mag'] != 0 else np.nan, axis=1
)

# Calculate summed_mean_div_metadata_contig_length, handling division by zero
df_ratio_analysis['summed_mean_div_metadata_contig_length'] = df_ratio_analysis.apply(
    lambda row: row['summed_mean_per_sample_mag'] / row['metadata_sum_contig_length']
                if row['metadata_sum_contig_length'] != 0 else np.nan, axis=1
)

print("\n5. All Ratio Calculations Complete.")
print(f"   Final analysis table created. Shape: {df_ratio_analysis.shape}")
print("   Head of final analysis table:\n", df_ratio_analysis.head())

print("\n   Descriptive statistics for the newly added ratios:")
print(df_ratio_analysis[['contig_length_ratio', 'summed_mean_div_sample_contig_length', 'summed_mean_div_metadata_contig_length']].describe())
print(f"   Number of NaN contig_length_ratios (due to zero metadata length): {df_ratio_analysis['contig_length_ratio'].isna().sum()}")


# --- 6. Save the new table ---
output_file_path = os.path.join(OUTPUT_DIR, 'extended_coverage_and_length_ratios_analysis_table.csv')
df_ratio_analysis.to_csv(output_file_path, index=False)
print(f"\nAnalysis table saved to: {output_file_path}")

print("\n--- Enhanced Ratio Analysis Script Complete ---")