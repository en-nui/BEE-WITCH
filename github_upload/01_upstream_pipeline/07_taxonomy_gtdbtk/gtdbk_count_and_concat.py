import os
import pandas as pd

# Define the base directory to search
base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/11_gtdbtk_output"

# Initialize counters and collections
gtdbtk_dir_count = 0
tsv_file_count = 0
dirs_without_tsv = 0
tsv_file_paths = []
processed_directories = set()

# Traverse the directory structure
for root, dirs, files in os.walk(base_dir):
    if "gtdbtk_output" in root:
        # Get the parent directory of the gtdbtk_output folder
        parent_directory = os.path.dirname(root)

        # Avoid processing the same parent directory multiple times
        if parent_directory in processed_directories:
            continue

        processed_directories.add(parent_directory)
        gtdbtk_dir_count += 1  # Count this gtdbtk_output directory

        # Prefer the classify subdirectory if it exists
        classify_tsv_path = os.path.join(root, "classify", "gtdbtk.bac120.summary.tsv")

        if os.path.isfile(classify_tsv_path):
            tsv_file_count += 1
            tsv_file_paths.append(classify_tsv_path)  # Use the file in the classify subdirectory
        else:
            dirs_without_tsv += 1

# Log the results
print(f"Total gtdbtk_output directories found: {gtdbtk_dir_count}")
print(f"Total gtdbtk.bac120.summary.tsv files found: {tsv_file_count}")
print(f"gtdbtk_output directories without tsv: {dirs_without_tsv}")

# Concatenate all summary.tsv files
output_file = "consolidated_summary.tsv"
df_list = []

for tsv_file in tsv_file_paths:
    # Read the tsv file
    try:
        df = pd.read_csv(tsv_file, sep="\t")

        # Extract parts of the directory structure
        parts = os.path.normpath(os.path.dirname(os.path.dirname(tsv_file))).split(os.sep)  # Go two levels up
        state, month, caste, colony = parts[-5], parts[-4], parts[-3], parts[-2]

        # Add extracted columns to the dataframe
        df["State"] = state
        df["Month"] = month
        df["Caste"] = caste
        df["Colony"] = colony

        df_list.append(df)
    except Exception as e:
        print(f"Error reading {tsv_file}: {e}")

# Combine all dataframes and deduplicate headers
if df_list:
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df.to_csv(output_file, sep="\t", index=False)
    print(f"Consolidated summary file written to {output_file}")
else:
    print("No TSV files were found to concatenate.")
