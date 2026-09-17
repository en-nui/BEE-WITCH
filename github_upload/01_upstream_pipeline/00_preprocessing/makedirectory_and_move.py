import os
import pandas as pd
import shutil

# Load the CSV file
csv_file = "Renamed_file_sheet.csv"
output_dir = "Organized_Files"  # Root output directory
os.makedirs(output_dir, exist_ok=True)  # Ensure the root directory exists

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_file)

# Iterate through each row of the DataFrame
for _, row in df.iterrows():
    renamed_file = row['renamed_file']
    state = str(row['State'])  # Ensure it's a string
    month = str(row['Month'])  # Ensure it's a string
    sample = str(row['Sample'])  # Ensure it's a string
    colony = str(row['Colony'])  # Ensure it's a string
    replicate = str(row['Replicate'])  # Ensure it's a string

    # Build the directory path
    replicate_dir = os.path.join(output_dir, state, month, sample, colony, replicate)
    os.makedirs(replicate_dir, exist_ok=True)  # Create the directories if they don't exist

    # Move the file into the appropriate directory
    source_path = os.path.join(renamed_file)  # Adjust if renamed_file contains a relative path
    destination_path = os.path.join(replicate_dir, os.path.basename(renamed_file))

    try:
        shutil.move(source_path, destination_path)
        print(f"Moved {renamed_file} to {replicate_dir}")
    except FileNotFoundError:
        print(f"File not found: {renamed_file}")
    except Exception as e:
        print(f"Error moving {renamed_file}: {e}")

print("Organizing and moving files completed!")

