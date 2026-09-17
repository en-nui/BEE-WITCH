import os
import pandas as pd
import re
from concurrent.futures import ProcessPoolExecutor

tmp_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/14_mMAG_dRep/tmp_directory"
metadata_file = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/14_mMAG_dRep/mMAG_contig_metadata_enriched.csv"
output_file = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/14_mMAG_dRep/organized_bacterial_bins.csv"


df_metadata = pd.read_csv(metadata_file)
metadata_dict = {
    (row['State'], row['Month'], row['Caste'], str(row['Colony']), row['contig_id']):
    (row['Family'], row['Genus'], row['Species'])
    for _, row in df_metadata.iterrows()
}

def parse_filename(filename):
    """Extracts metadata from filename."""
    match = re.match(r"(\w+)_(\w+)_(\w+)_(\d+)_final_bins_bin_(\d+)\.fna", filename)
    if match:
        state, month, caste, colony, bin_id = match.groups()
        return state, month, caste, colony, f"bin_{bin_id}", filename
    return None

def process_fna_file(file):
    """Processes a single .fna file and extracts relevant data."""
    file_path = os.path.join(tmp_dir, file)
    parsed = parse_filename(file)
    if not parsed:
        return []
    
    state, month, caste, colony, bin_id, source_file = parsed
    data = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_id = line.strip().lstrip('>')
                family, genus, species = metadata_dict.get((state, month, caste, colony, contig_id), ("Unknown", "Unknown", "Unknown"))
                data.append([contig_id, state, month, caste, colony, bin_id, source_file, family, genus, species])
    
    return data

# Use multiprocessing to speed up file processing
data = []
fna_files = [f for f in os.listdir(tmp_dir) if f.endswith(".fna")]

with ProcessPoolExecutor() as executor:
    results = executor.map(process_fna_file, fna_files)
    for result in results:
        data.extend(result)

# Create DataFrame and save
columns = ["contig_id", "State", "Month", "Caste", "Colony", "bin_id", "source_file", "Family", "Genus", "Species"]
df_output = pd.DataFrame(data, columns=columns)
df_output.to_csv(output_file, index=False)

print(f"CSV saved to {output_file}")

