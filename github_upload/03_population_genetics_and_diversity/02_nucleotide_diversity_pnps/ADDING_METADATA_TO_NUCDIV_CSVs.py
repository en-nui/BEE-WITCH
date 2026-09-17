import pandas as pd
import os
import numpy as np # Import numpy for np.nan

# Define file paths
updated_annotation_path = '/home/robinch/projects/BEE-WITCH/22_anvio/UPDATED_ANNOTATION_with_full_metadata_20250416.csv'
nucdiv_gene_filtered_path = '/home/robinch/projects/BEE-WITCH/22_anvio/all_nucdiv_by_gene_RETAINED_filtered.csv'

# Load the updated annotation data and prepare it for merging
print("Loading UPDATED_ANNOTATION_with_full_metadata_20250416.csv...")
try:
    updated_annotation_df = pd.read_csv(updated_annotation_path)
    print("Loaded UPDATED_ANNOTATION_with_full_metadata_20250416.csv successfully.")

    # Select relevant columns for functional category
    required_cols_annotation = ['updated_contig_id', 'gene_position', 'functional_category']
    if not all(col in updated_annotation_df.columns for col in required_cols_annotation):
        missing = [col for col in required_cols_annotation if col not in updated_annotation_df.columns]
        print(f"Error: Missing required columns in {updated_annotation_path} for functional category merging: {missing}")
        exit()

    annotation_info_df = updated_annotation_df[required_cols_annotation].copy()

    # Ensure updated_contig_id is string and strip whitespace
    print("Cleaning 'updated_contig_id' in annotation data...")
    annotation_info_df['updated_contig_id'] = annotation_info_df['updated_contig_id'].astype(str).str.strip()

    # Ensure gene_position is integer and handle errors
    print("Converting 'gene_position' to integer in annotation data...")
    annotation_info_df['gene_position'] = pd.to_numeric(annotation_info_df['gene_position'], errors='coerce')
    # Drop rows where gene_position could not be converted
    initial_annotation_rows = len(annotation_info_df)
    annotation_info_df.dropna(subset=['gene_position'], inplace=True)
    rows_removed_annotation = initial_annotation_rows - len(annotation_info_df)
    if rows_removed_annotation > 0:
        print(f"Warning: Removed {rows_removed_annotation} rows from annotation data due to non-numeric 'gene_position'.")

    # Convert to integer type after dropping NaNs
    annotation_info_df['gene_position'] = annotation_info_df['gene_position'].astype(int)


except FileNotFoundError:
    print(f"Error: Annotation file not found at {updated_annotation_path}")
    exit()
except Exception as e:
    print(f"An error occurred while processing the annotation file: {e}")
    exit()


def add_annotation_to_nucdiv(nucdiv_file_path, annotation_data_df):
    """
    Loads a filtered nucleotide diversity file, merges with annotation data, and saves.
    """
    print(f"\nProcessing file: {os.path.basename(nucdiv_file_path)}")
    try:
        nucdiv_df = pd.read_csv(nucdiv_file_path)
        initial_cols = nucdiv_df.columns.tolist()
        print(f"Initial columns in {os.path.basename(nucdiv_file_path)}: {initial_cols}")

        # Ensure required columns for merging exist in nucdiv_df
        required_cols_nucdiv = ['updated_contig_id', 'gene_position'] # Use gene_position as per your renaming
        if not all(col in nucdiv_df.columns for col in required_cols_nucdiv):
            missing = [col for col in required_cols_nucdiv if col not in nucdiv_df.columns]
            print(f"Error: Missing required merge columns in {os.path.basename(nucdiv_file_path)}: {missing}. Skipping.")
            return

        # Ensure updated_contig_id is string and strip whitespace
        print(f"Cleaning 'updated_contig_id' in {os.path.basename(nucdiv_file_path)}...")
        nucdiv_df['updated_contig_id'] = nucdiv_df['updated_contig_id'].astype(str).str.strip()

        # Ensure gene_position is integer and handle errors
        print(f"Converting 'gene_position' to integer in {os.path.basename(nucdiv_file_path)}...")
        nucdiv_df['gene_position'] = pd.to_numeric(nucdiv_df['gene_position'], errors='coerce')
         # Drop rows where gene_position could not be converted
        initial_nucdiv_rows_for_merge = len(nucdiv_df)
        nucdiv_df.dropna(subset=['gene_position'], inplace=True)
        rows_removed_nucdiv = initial_nucdiv_rows_for_merge - len(nucdiv_df)
        if rows_removed_nucdiv > 0:
            print(f"Warning: Removed {rows_removed_nucdiv} rows from {os.path.basename(nucdiv_file_path)} due to non-numeric 'gene_position' before merge.")

        # Convert to integer type after dropping NaNs
        nucdiv_df['gene_position'] = nucdiv_df['gene_position'].astype(int)


        # --- Merge with annotation data ---
        print("Merging with annotation data...")

        # Check number of rows before merge (after gene_position cleaning)
        print(f"Rows in {os.path.basename(nucdiv_file_path)} ready for merge: {len(nucdiv_df)}")
        print(f"Rows in annotation_data_df ready for merge: {len(annotation_data_df)}")


        # Perform a left merge to add annotation information
        updated_nucdiv_df = pd.merge(nucdiv_df,
                                     annotation_data_df,
                                     on=['updated_contig_id', 'gene_position'], # Merge on both columns directly
                                     how='left')


        # Check how many rows got a functional_category after the merge
        rows_with_functional_category = updated_nucdiv_df['functional_category'].notna().sum()
        print(f"Rows in {os.path.basename(nucdiv_file_path)} after merge with functional_category: {rows_with_functional_category}")


        final_cols = updated_nucdiv_df.columns.tolist()
        print(f"Final columns in {os.path.basename(nucdiv_file_path)} after merge: {final_cols}")

        # Save the updated DataFrame
        output_file_path = nucdiv_file_path.replace('_RETAINED_filtered.csv', '_RETAINED_annotated.csv') # Save to a new file
        updated_nucdiv_df.to_csv(output_file_path, index=False)
        print(f"Updated data saved to: {output_file_path}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {nucdiv_file_path}. Skipping.")
    except Exception as e:
        print(f"An error occurred while processing {os.path.basename(nucdiv_file_path)}: {e}. Skipping.")


# Process only the gene-level filtered nucleotide diversity file
add_annotation_to_nucdiv(nucdiv_gene_filtered_path, annotation_info_df)

print("\nFinished adding annotation data to the gene-level filtered nucleotide diversity file.")
