import os
import logging
import subprocess
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Define base directories
base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH"
trimmed_reads_dir = os.path.join(base_dir, "01_TrimmedReads")
assemblies_dir = os.path.join(base_dir, "02a_LarvalCoassembly")

def find_larvae_directories(base_dir):
    """
    Identifies trimmed reads for larval samples in each state and concatenates them for coassembly.
    Returns a dictionary of state -> (combined_R1_path, combined_R2_path).
    """
    larvae_dirs = {}
    for state_dir in Path(base_dir).iterdir():
        if not state_dir.is_dir():
            continue
        state = state_dir.name
        logging.info(f"Processing state directory: {state_dir}")

        combined_r1_path = os.path.join(state_dir, f"{state}_combined_R1.fastq.gz")
        combined_r2_path = os.path.join(state_dir, f"{state}_combined_R2.fastq.gz")

        # Traverse through the hierarchy to find reads
        with open(combined_r1_path, "wb") as r1_out, open(combined_r2_path, "wb") as r2_out:
            for month_dir in state_dir.iterdir():
                larvae_dir = month_dir / "Larvae"
                if not larvae_dir.is_dir():
                    continue
                logging.info(f"Looking for reads in: {larvae_dir}")

                for colony_dir in larvae_dir.iterdir():
                    if not colony_dir.is_dir():
                        continue
                    for sub_dir in colony_dir.iterdir():
                        if not sub_dir.is_dir():
                            continue
                        # Find R1 and R2 reads
                        r1_files = list(sub_dir.glob("*_R1.trimmed.fastq.gz"))
                        r2_files = list(sub_dir.glob("*_R2.trimmed.fastq.gz"))

                        if len(r1_files) != len(r2_files):
                            logging.warning(f"Mismatch in R1/R2 files in {sub_dir}")
                            continue

                        # Concatenate R1 and R2 reads
                        for r1_file, r2_file in zip(r1_files, r2_files):
                            logging.info(f"Adding {r1_file} to {combined_r1_path}")
                            with open(r1_file, "rb") as r1_in:
                                r1_out.write(r1_in.read())
                            logging.info(f"Adding {r2_file} to {combined_r2_path}")
                            with open(r2_file, "rb") as r2_in:
                                r2_out.write(r2_in.read())

        larvae_dirs[state] = (combined_r1_path, combined_r2_path)
        logging.info(f"Finished combining reads for state: {state}")
    return larvae_dirs

def run_megahit(state, r1_path, r2_path):
    """
    Runs MEGAHIT on concatenated R1 and R2 reads for a given state.
    """
    output_dir = os.path.join(assemblies_dir, state)
    #os.makedirs(output_dir, exist_ok=True)

    megahit_cmd = [
        "megahit",
        "-1", r1_path,
        "-2", r2_path,
        "-o", output_dir,
        "--presets", "meta-sensitive",
        "-t", "24",
        "--continue",
    ]

    try:
        subprocess.run(megahit_cmd, check=True)
        logging.info(f"MEGAHIT completed for state: {state}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running MEGAHIT for state {state}: {e}")

def process_larval_coassemblies():
    """
    Main function to process larval coassemblies.
    """
    larvae_dirs = find_larvae_directories(trimmed_reads_dir)
    if not larvae_dirs:
        logging.warning("No larval reads found for coassembly.")
        return

    logging.info(f"Starting coassembly for {len(larvae_dirs)} states.")
    for state, (r1_path, r2_path) in larvae_dirs.items():
        logging.info(f"Starting coassembly for state: {state}")
        run_megahit(state, r1_path, r2_path)

    logging.info("All coassemblies completed.")

if __name__ == "__main__":
    process_larval_coassemblies()

