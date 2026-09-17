import os
import logging
import datetime
from concurrent.futures import ThreadPoolExecutor, wait
import pandas as pd
from Bio import SeqIO

# 1. Define paths (UPDATED)
output_base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/16_mMAG_annotation"
master_fasta_file = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/15_MAG_dRep/renamed_concat_genomes/concat_mMAG.fna"
contig_info_csv = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_HOME_NO_TOUCH_20250211.csv"

# 2. Set up logging (No change needed)
log_file = os.path.join(output_base_dir, "annotation_process.log")
logging.basicConfig(filename=log_file, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Starting colony-specific mMAG contig extraction and parallel DRAM annotation process.")

# 3. Load Contig Info CSV (NEW)
try:
    contig_df = pd.read_csv(contig_info_csv)
    logging.info(f"Contig information CSV loaded successfully from: {contig_info_csv}")
except FileNotFoundError:
    logging.error(f"Error: Contig information CSV file not found at: {contig_info_csv}. Exiting.")
    print(f"Error: Contig information CSV file not found at: {contig_info_csv}. Check log file for details.")
    exit()

# 4. Load Master FASTA File into a dictionary for fast lookup (NEW)
contig_sequences = {}
try:
    with open(master_fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig_sequences[record.id] = record.seq
    logging.info(f"Master FASTA file loaded into memory: {master_fasta_file}")
    logging.info(f"Number of contigs in master FASTA: {len(contig_sequences)}")
except FileNotFoundError:
    logging.error(f"Error: Master FASTA file not found at: {master_fasta_file}. Exiting.")
    print(f"Error: Master FASTA file not found at: {master_fasta_file}. Check log file for details.")
    exit()

# 5. Get Unique Colonies and Process Each (UPDATED)
unique_colonies = contig_df['Colony'].unique()
dram_commands = []

for colony_id in unique_colonies:
    colony_folder = os.path.join(output_base_dir, f"colony_{colony_id}")
    concat_filepath = os.path.join(colony_folder, "tmp_concat.fna")
    dram_output_dir = os.path.join(colony_folder, "DRAM_output")
    os.makedirs(colony_folder, exist_ok=True)
    #os.makedirs(dram_output_dir, exist_ok=True) # Created when DRAM is executed if needed

    logging.info(f"Starting contig extraction for colony: {colony_id}")

    # Filter DataFrame for current colony and mMAG type
    colony_contigs_df = contig_df[(contig_df['Colony'] == colony_id) & (contig_df['type'] == 'mMAG')]
    colony_contig_names_csv = colony_contigs_df['updated_contig_id'].tolist() # Contigs listed in CSV
    num_contigs_csv = len(colony_contig_names_csv)
    logging.info(f"Found {num_contigs_csv} mMAG contigs for colony {colony_id} in CSV.")


    concat_start_time = datetime.datetime.now()
    contigs_copied_count = 0 # Counter for actually copied contigs
    with open(concat_filepath, 'w') as outfile:
        for contig_name in colony_contig_names_csv:
            if contig_name in contig_sequences:
                sequence = contig_sequences[contig_name]
                outfile.write(f">{contig_name}\n{sequence}\n")
                contigs_copied_count += 1
            else:
                logging.warning(f"Contig '{contig_name}' from CSV not found in master FASTA file!")

    concat_end_time = datetime.datetime.now()
    concat_duration = concat_end_time - concat_start_time
    logging.info(f"Finished contig extraction for colony: {colony_id}. Duration: {concat_duration}")
    logging.info(f"Copied {contigs_copied_count} contigs to {concat_filepath}")
    print(f"Colony {colony_id}: Found {num_contigs_csv} mMAG contigs in CSV, Copied {contigs_copied_count} contigs to {concat_filepath}")


    # Construct DRAM command (No change needed in command itself)
    dram_command = f"DRAM.py annotate -i {concat_filepath} -o {dram_output_dir}"
    logging.info(f"DRAM command for colony {colony_id}: {dram_command}")
    print(f"DRAM command for colony {colony_id}: {dram_command}")
    dram_commands.append(dram_command)

# 6. Parallel execution of DRAM commands (No change needed)
def execute_dram_command(command):
    logging.info(f"Thread {os.getpid()}: Executing DRAM command: {command}")
    start_time = datetime.datetime.now()
    os.system(command) # Or subprocess.run for better error handling if needed
    end_time = datetime.datetime.now()
    duration = end_time - start_time
    logging.info(f"Thread {os.getpid()}: Finished DRAM command. Duration: {duration}")

max_workers = 4
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(execute_dram_command, cmd) for cmd in dram_commands]
    wait(futures) # Wait for all DRAM commands to complete

logging.info("Colony-specific mMAG contig extraction and parallel DRAM annotation process finished.")
print("Colony-specific mMAG contig extraction and parallel DRAM annotation process finished. Check annotation_process.log for details and DRAM commands.")
