#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 22 14:12:44 2025

@author: robinch
"""

import pandas as pd
import os


df = pd.read_csv('/home/robinch/projects/BEE-WITCH/16_rMAG_mapping/colony_rMAG_mapping_with_metadata/concat_input/concat_annotated_summaries_for_R.csv')
df_mMAG = df[df['type']=="mMAG"]

# Directory to save the individual colony DataFrames
output_dir = '/home/robinch/projects/BEE-WITCH/29_mMAG_analyses/colony_data_splits'
os.makedirs(output_dir, exist_ok=True) # Create the directory if it doesn't exist

# --- 1. Generate Summary Statistics ---
print("\n--- Summary Statistics ---")

# Number of unique MAG_id row values
num_unique_mag_ids = df_mMAG['MAG_id'].nunique()
print(f"Number of unique MAG_id values: {num_unique_mag_ids}")

# Number of unique Species row values
num_unique_species = df_mMAG['Species'].nunique()
print(f"Number of unique Species values: {num_unique_species}")

# Number of unique Genus values
num_unique_genus = df_mMAG['Genus'].nunique()
print(f"Number of unique Genus values: {num_unique_genus}")

# Within each unique Genus, the number of unique MAG_id values
print("\nNumber of unique MAG_id values within each Genus:")
mag_ids_per_genus = df_mMAG.groupby('Genus')['MAG_id'].nunique()
print(mag_ids_per_genus.to_string()) # .to_string() for better display of all rows

# Within each unique Species, the number of unique MAG_id values
print("\nNumber of unique MAG_id values within each Species:")
mag_ids_per_species = df_mMAG.groupby('Species')['MAG_id'].nunique()
print(mag_ids_per_species.to_string()) # .to_string() for better display of all rows


# Within each unique primary_cluster, the number of unique primary_cluster values
print("\nNumber of unique MAG_id values within each Species:")
primary_cluster_per_species = df_mMAG.groupby('Species')['primary_cluster'].nunique()
print(primary_cluster_per_species.to_string()) # .to_string() for better display of all rows

# --- 2. Split Data by Colony (Memory Efficiently) ---
print("\n--- Splitting Data by Colony ---")

unique_colonies = df_mMAG['Colony'].unique()
print(f"Found {len(unique_colonies)} unique Colony values.")

for colony_val in unique_colonies:
    # Filter the main DataFrame for the current colony.
    # This creates a *view* or a *copy* depending on pandas internal logic.
    # For saving to CSV, a copy is implicitly made when writing, so it's fine.
    colony_df = df_mMAG[df_mMAG['Colony'] == colony_val]

    # Sanitize colony_val for filename (replace characters that might cause issues)
    safe_colony_val = str(colony_val).replace('/', '_').replace('\\', '_').replace(' ', '_')
    output_filepath = os.path.join(output_dir, f'colony_{safe_colony_val}.csv')

    # Save the filtered DataFrame to a new CSV file
    colony_df.to_csv(output_filepath, index=False)
    print(f"Saved data for {colony_val} to {output_filepath}")

print(f"\nAll colony data has been split and saved to the '{output_dir}' directory.")
print("Memory efficiency achieved by saving to disk rather than holding all sub-DataFrames in memory.")
