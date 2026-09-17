import pandas as pd
import os

# Define file paths
retained_mags_path = '/home/robinch/projects/BEE-WITCH/21_final_r_vMAG_set/all_colonies_RETAINED_r_vMAGs.csv'
nucdiv_gene_path = '/home/robinch/projects/BEE-WITCH/22_anvio/all_nucdiv_by_gene.csv'
nucdiv_mag_path = '/home/robinch/projects/BEE-WITCH/22_anvio/all_nucdiv_by_mag.csv'

# Load the retained MAGs data
try:
    retained_mags_df = pd.read_csv(retained_mags_path)
    print("Loaded all_colonies_RETAINED_r_vMAGs.csv successfully.")
    # Ensure relevant columns are strings for merging, if necessary
    for col in ['MAG_id', 'State', 'Month', 'Caste', 'Colony', 'Rep']:
         if col in retained_mags_df.columns:
             retained_mags_df[col] = retained_mags_df[col].astype(str)
         else:
             print(f"Warning: Column '{col}' not found in {retained_mags_path}")

except FileNotFoundError:
    print(f"Error: Retained MAGs file not found at {retained_mags_path}")
    exit()

def filter_nucdiv_file(file_path, retained_df):
    """
    Loads a nucleotide diversity file, splits sample_id, filters based on retained data,
    reports removed rows, and saves the filtered data.
    """
    print(f"\nProcessing file: {file_path}")
    try:
        nucdiv_df = pd.read_csv(file_path)
        initial_rows = len(nucdiv_df)
        print(f"Initial rows in {os.path.basename(file_path)}: {initial_rows}")

        # Check if 'sample_id' and 'MAG_id' columns exist
        if 'sample_id' not in nucdiv_df.columns:
            print(f"Error: 'sample_id' column not found in {os.path.basename(file_path)}. Skipping.")
            return
        if 'MAG_id' not in nucdiv_df.columns:
             print(f"Error: 'MAG_id' column not found in {os.path.basename(file_path)}. Skipping.")
             return


        # Extract metadata from sample_id
        print("Extracting metadata from 'sample_id'...")
        try:
            # Split sample_id and create new columns
            nucdiv_df[['State', 'Month', 'Caste', 'Colony', 'Rep']] = nucdiv_df['sample_id'].str.split('_', expand=True)

            # Ensure extracted columns are strings for merging
            for col in ['State', 'Month', 'Caste', 'Colony', 'Rep']:
                nucdiv_df[col] = nucdiv_df[col].astype(str)


        except Exception as e:
            print(f"Error splitting 'sample_id' in {os.path.basename(file_path)}: {e}. Skipping.")
            return

        # Ensure MAG_id is string for merging
        nucdiv_df['MAG_id'] = nucdiv_df['MAG_id'].astype(str)

        # --- Filter based on retained data ---
        print("Filtering based on retained MAGs data...")
        # Define the columns to merge on
        merge_cols = ['MAG_id', 'State', 'Month', 'Caste', 'Colony', 'Rep']

        # Ensure all merge columns exist in both dataframes
        if not all(col in nucdiv_df.columns for col in merge_cols):
            missing_cols = [col for col in merge_cols if col not in nucdiv_df.columns]
            print(f"Error: Missing merge columns in {os.path.basename(file_path)}: {missing_cols}. Skipping.")
            return
        if not all(col in retained_df.columns for col in merge_cols):
             missing_cols = [col for col in merge_cols if col not in retained_df.columns]
             print(f"Error: Missing merge columns in {os.path.basename(retained_mags_path)}: {missing_cols}. Skipping.")
             return


        # Perform inner merge to keep only matching rows
        filtered_nucdiv_df = pd.merge(nucdiv_df,
                                     retained_df[merge_cols], # Select only the merge columns from retained_df
                                     on=merge_cols,
                                     how='inner')

        final_rows = len(filtered_nucdiv_df)
        rows_removed = initial_rows - final_rows
        print(f"Filtered rows in {os.path.basename(file_path)}: {final_rows}")
        print(f"Number of rows removed during filtering: {rows_removed}")

        # Save the filtered DataFrame
        output_file_path = file_path.replace('.csv', '_RETAINED_filtered.csv') # Save to a new file
        filtered_nucdiv_df.to_csv(output_file_path, index=False)
        print(f"Filtered data saved to: {output_file_path}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {file_path}. Skipping.")
    except Exception as e:
        print(f"An error occurred while processing {os.path.basename(file_path)}: {e}. Skipping.")


# Process the two nucleotide diversity files
filter_nucdiv_file(nucdiv_gene_path, retained_mags_df)
filter_nucdiv_file(nucdiv_mag_path, retained_mags_df)

print("\nFinished filtering nucleotide diversity files.")
