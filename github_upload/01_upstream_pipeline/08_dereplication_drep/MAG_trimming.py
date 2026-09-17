import csv
from Bio import SeqIO

combined_contig_metadata_path = 'RENAMED_CONTIG_HOME_NO_TOUCH_20250211.csv' # Assuming this is in the same directory
concat_vmag_path = '/home/robinch/projects/BEE-WITCH/15_MAG_dRep/renamed_concat_genomes/concat_vMAG.fna'
concat_mmag_path = '/home/robinch/projects/BEE-WITCH/15_MAG_dRep/renamed_concat_genomes/concat_mMAG.fna'
output_trimmed_vmag_path = 'trimmed_concat_vMAG.fna'
output_trimmed_mmag_path = 'trimmed_concat_mMAG.fna'
log_file_path = 'mag_trimming_log.log'


# 1. Load combined contig metadata and identify representative contig IDs and counts
representative_contig_ids = set()
unique_r_mag_ids = set()
total_representative_contigs_metadata = 0

try:
    with open(combined_contig_metadata_path, 'r', newline='') as metadata_infile:
        metadata_reader = csv.DictReader(metadata_infile)
        for row in metadata_reader:
            if row['MAG_id'] == row['r_MAG_id']: # Check if MAG_id is representative
                representative_contig_ids.add(row['updated_contig_id'])
                unique_r_mag_ids.add(row['r_MAG_id'])
                total_representative_contigs_metadata += 1 # Count contigs associated with representative MAGs
except FileNotFoundError:
    print(f"Error: Combined contig metadata file '{combined_contig_metadata_path}' not found.")
    exit()

# 2. Trim concat_vMAG.fna
trimmed_vmag_records = []
vmag_contigs_kept = 0
vmag_contigs_discarded = 0

try:
    with open(concat_vmag_path, 'r') as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in representative_contig_ids:
                trimmed_vmag_records.append(record)
                vmag_contigs_kept += 1
            else:
                vmag_contigs_discarded += 1
except FileNotFoundError:
    print(f"Warning: vMAG concatenated file '{concat_vmag_path}' not found. Skipping vMAG trimming.")

# 3. Trim concat_mMAG.fna
trimmed_mmag_records = []
mmag_contigs_kept = 0
mmag_contigs_discarded = 0

try:
    with open(concat_mmag_path, 'r') as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in representative_contig_ids:
                trimmed_mmag_records.append(record)
                mmag_contigs_kept += 1
            else:
                mmag_contigs_discarded += 1
except FileNotFoundError:
    print(f"Warning: mMAG concatenated file '{concat_mmag_path}' not found. Skipping mMAG trimming.")


# 4. Write trimmed FASTA files and log
with open(log_file_path, 'w') as logfile:
    logfile.write("Trimming Summary:\n\n")
    logfile.write(f"Unique representative r_MAG_IDs found in metadata: {len(unique_r_mag_ids)}\n")
    logfile.write(f"Total contigs associated with representative r_MAG_IDs (in metadata): {total_representative_contigs_metadata}\n\n")


    if trimmed_vmag_records:
        SeqIO.write(trimmed_vmag_records, output_trimmed_vmag_path, "fasta")
        logfile.write(f"Trimmed vMAG contigs saved to '{output_trimmed_vmag_path}'. Kept: {vmag_contigs_kept}, Discarded: {vmag_contigs_discarded}, Headers in trimmed file: {len(trimmed_vmag_records)}\n")
        print(f"Trimmed vMAG contigs saved to '{output_trimmed_vmag_path}'")
    else:
        logfile.write("No representative vMAG contigs found or vMAG input file missing. Output vMAG file not created.\n")
        print("No representative vMAG contigs found or vMAG input file missing. Output vMAG file not created.")

    if trimmed_mmag_records:
        SeqIO.write(trimmed_mmag_records, output_trimmed_mmag_path, "fasta")
        logfile.write(f"Trimmed mMAG contigs saved to '{output_trimmed_mmag_path}'. Kept: {mmag_contigs_kept}, Discarded: {mmag_contigs_discarded}, Headers in trimmed file: {len(trimmed_mmag_records)}\n")
        print(f"Trimmed mMAG contigs saved to '{output_trimmed_mmag_path}'")
    else:
        logfile.write("No representative mMAG contigs found or mMAG input file missing. Output mMAG file not created.\n")
        print("No representative mMAG contigs found or mMAG input file missing. Output mMAG file not created.")

print(f"Trimming summary logged to '{log_file_path}'")
