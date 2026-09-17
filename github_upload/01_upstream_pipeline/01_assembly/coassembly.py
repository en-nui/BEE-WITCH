import os
import subprocess
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Define input and output base directories
trimmed_reads_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/01_TrimmedReads"
assemblies_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/03_Assemblies"

def cleanup_temp_files(base_dir):
    """
    Removes temp_R1.fq.gz and temp_R2.fq.gz files from the trimmed reads directory.
    """
    logging.info("Starting cleanup of temporary files.")
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file in {"temp_R1.fq.gz", "temp_R2.fq.gz"}:
                temp_file_path = os.path.join(root, file)
                try:
                    os.remove(temp_file_path)
                    logging.info(f"Removed temporary file: {temp_file_path}")
                except Exception as e:
                    logging.error(f"Error removing file {temp_file_path}: {e}")
    logging.info("Cleanup of temporary files completed.")

def find_coassembly_directories(base_dir):
    """
    Traverses directory structure to find parent directories containing trimmed fastq.gz files.
    """
    coassembly_dirs = {}
    
    for root, dirs, files in os.walk(base_dir):
        # Filter R1/R2 trimmed fastq files
        r1_files = sorted(f for f in files if f.endswith("_R1.trimmed.fastq.gz"))
        r2_files = sorted(f for f in files if f.endswith("_R2.trimmed.fastq.gz"))
        
        # If R1 and R2 files are found, add their parent directory to the coassembly list
        if r1_files and r2_files:
            parent_dir = Path(root).parent
            if parent_dir not in coassembly_dirs:
                coassembly_dirs[parent_dir] = []
            coassembly_dirs[parent_dir].extend(
                [(os.path.join(root, r1), os.path.join(root, r2)) for r1, r2 in zip(r1_files, r2_files)]
            )

    logging.info(f"Found {len(coassembly_dirs)} directories with R1/R2 pairs for coassembly.")
    return coassembly_dirs

def run_megahit_with_concat(input_dir, assembly_dir, file_pairs):
    """
    Runs MEGAHIT coassembly on R1/R2 pairs within a given directory.
    """
    logging.info(f"Running coassembly for: {input_dir}")

    # Ensure the output directory exists
    os.makedirs(assembly_dir, exist_ok=True)

    # Temporary files for concatenated inputs
    temp_r1 = os.path.join(input_dir, "temp_R1.fq.gz")
    temp_r2 = os.path.join(input_dir, "temp_R2.fq.gz")

    # Concatenate all R1 and R2 files
    with open(temp_r1, "wb") as r1_out, open(temp_r2, "wb") as r2_out:
        for r1, r2 in file_pairs:
            logging.info(f"Adding {r1} to {temp_r1}")
            with open(r1, "rb") as r1_in:
                r1_out.write(r1_in.read())
            logging.info(f"Adding {r2} to {temp_r2}")
            with open(r2, "rb") as r2_in:
                r2_out.write(r2_in.read())

    # Run MEGAHIT
    megahit_cmd = [
        "megahit",
        "-1", temp_r1,
        "-2", temp_r2,
        "-o", assembly_dir,
        "--presets", "meta-sensitive",
    ]
    try:
        subprocess.run(megahit_cmd, check=True)
        logging.info(f"Coassembly completed for: {input_dir}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during MEGAHIT coassembly for {input_dir}: {e}")

def process_coassemblies():
    """
    Main function to process coassemblies in parallel.
    """
    # Cleanup temporary files first
    cleanup_temp_files(trimmed_reads_dir)
    
    # Find directories with R1/R2 pairs for coassembly
    coassembly_dirs = find_coassembly_directories(trimmed_reads_dir)
    
    # Prepare input and output directories
    tasks = []
    for parent_dir, file_pairs in coassembly_dirs.items():
        assembly_dir = os.path.join(assemblies_dir, os.path.relpath(parent_dir, trimmed_reads_dir))
        tasks.append((parent_dir, assembly_dir, file_pairs))
    
    # Run MEGAHIT coassemblies in parallel
    with ProcessPoolExecutor(max_workers=24) as executor:
        for parent_dir, assembly_dir, file_pairs in tasks:
            executor.submit(run_megahit_with_concat, parent_dir, assembly_dir, file_pairs)

    logging.info("All coassemblies completed.")

if __name__ == "__main__":
    process_coassemblies()

