import os
import pandas as pd

# Load the CSV file
csv_file = "Renamed_file_sheet.csv"
df = pd.read_csv(csv_file)

# Ensure the necessary columns are present
if not {"original_file_name", "renamed_file"}.issubset(df.columns):
    raise ValueError("CSV file must contain 'original_file_name' and 'renamed_file' columns.")

# Directory containing the files to rename
input_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/ILMN_2356_CGB_IUB_GSF3924_Nov2024"

# Rename files
for _, row in df.iterrows():
    original_path = os.path.join(input_dir, row["original_file_name"])
    renamed_path = os.path.join(input_dir, row["renamed_file"])

    # Check if the original file exists
    if os.path.exists(original_path):
        os.rename(original_path, renamed_path)
        print(f"Renamed: {original_path} -> {renamed_path}")
    else:
        print(f"File not found: {original_path}")
