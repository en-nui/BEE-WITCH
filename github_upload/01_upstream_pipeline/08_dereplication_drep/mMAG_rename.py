import csv
import os
import re
from Bio import SeqIO

mmag_fna_directory = '/home/robinch/projects/BEE-WITCH/14_mMAG_dRep/tmp_directory/'
combined_contig_metadata_path = 'RENAMED_CONTIG_HOME_NO_TOUCH_20250211.csv' # Assuming this is in the same directory
output_concat_mmag_path = 'concat_mMAG.fna'
log_file_path = 'mmag_rename_log.log'

def rename_mmag_contigs(mmag_fna_file, metadata_for_mmag, logfile):
    """
    Renames contig headers in a mMAG FNA file based on metadata.
    Logs any mismatches or errors.
    Returns a list of SeqIO records with renamed headers.
    """
    renamed_records = []
    contig_id_map = {row['original_contig_id']: row['updated_contig_id']
                      for row in metadata_for_mmag if row['original_contig_id'] != 'No contig_header'} # Create quick lookup
    contigs_in_fna = 0

    try:
        for record in SeqIO.parse(mmag_fna_file, "fasta"):
            contigs_in_fna += 1
            original_header = str(record.id) # Ensure it's a string
            updated_header = contig_id_map.get(original_header)

            if updated_header:
                record.id = updated_header
                record.description = "" # Clear description
                renamed_records.append(record)
            else:
                logfile.write(f"Warning: Contig header '{original_header}' from '{mmag_fna_file}' not found in metadata.\n")

    except Exception as e:
        logfile.write(f"Error processing FNA file '{mmag_fna_file}': {e}\n")
        return None # Indicate error

    metadata_contig_count = len(metadata_for_mmag)
    if contigs_in_fna != metadata_contig_count:
        logfile.write(f"Warning: Contig count mismatch for '{mmag_fna_file}'. FNA file has {contigs_in_fna} contigs, metadata expects {metadata_contig_count}.\n")

    return renamed_records

def parse_filename_metadata(filename):
    """Parses metadata from the filename. Returns a dictionary or None if parsing fails."""
    pattern = r"([^_]+)_([^_]+)_([^_]+)_([^_]+)_final_bins_(.*)\.fna"
    match = re.match(pattern, filename)
    if match:
        return {
            'State': match.group(1),
            'Month': match.group(2),
            'Caste': match.group(3),
            'Colony': match.group(4),
            'bin_id': match.group(5)
        }
    else:
        return None

def construct_filename(state, month, caste, colony, user_genome):
    """Constructs the filename from metadata components."""
    return f"{state}_{month}_{caste}_{colony}_final_bins_{user_genome}.fna"

# 1. Load combined contig metadata and filter for mMAGs
mmag_metadata = {} # Dictionary to store mMAG metadata, keyed by filename metadata components
try:
    with open(combined_contig_metadata_path, 'r', newline='') as metadata_infile:
        metadata_reader = csv.DictReader(metadata_infile)
        for row in metadata_reader:
            if row['type'] == 'mMAG':
                metadata_key = (row['State'], row['Month'], row['Caste'], row['Colony'], row['bin_id'])
                if metadata_key not in mmag_metadata:
                    mmag_metadata[metadata_key] = []
                mmag_metadata[metadata_key].append(row)
except FileNotFoundError:
    print(f"Error: Combined contig metadata file '{combined_contig_metadata_path}' not found.")
    exit()

# 2. Process mMAG FNA files and rename contigs
renamed_mmag_records = []
mmag_files_processed = 0
mmag_files_not_found = 0


with open(log_file_path, 'w') as logfile: # Open log file for writing
    for filename in os.listdir(mmag_fna_directory):
        if filename.endswith(".fna"): # Process all .fna files in the mMAG directory
            mmag_files_processed += 1
            mmag_fna_file_path = os.path.join(mmag_fna_directory, filename)
            filename_metadata = parse_filename_metadata(filename)

            if filename_metadata:
                metadata_key = (filename_metadata['State'], filename_metadata['Month'], filename_metadata['Caste'], filename_metadata['Colony'], filename_metadata['bin_id'])
                metadata_for_current_mmag = mmag_metadata.get(metadata_key)

                if metadata_for_current_mmag:
                    renamed_records = rename_mmag_contigs(mmag_fna_file_path, metadata_for_current_mmag, logfile)
                    if renamed_records:
                        renamed_mmag_records.extend(renamed_records)
                else:
                    logfile.write(f"Warning: No metadata found in combined sheet for filename metadata: {filename_metadata} (filename: '{filename}').\n")

            else:
                logfile.write(f"Warning: Could not parse filename metadata for '{filename}'. Skipping FNA file.\n")


# 3. Concatenate all renamed mMAG records
if renamed_mmag_records:
    SeqIO.write(renamed_mmag_records, output_concat_mmag_path, "fasta")
    print(f"Successfully renamed and concatenated mMAG contigs to '{output_concat_mmag_path}'")
else:
    print("No mMAG contigs were processed or renamed. Output file not created.")

print(f"mMAG files processed: {mmag_files_processed}")
print(f"Check '{log_file_path}' for warnings and errors.")
