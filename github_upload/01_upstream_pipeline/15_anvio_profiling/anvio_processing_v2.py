import os
import glob
import logging
import subprocess
import argparse
import multiprocessing

# Define default values (can be overridden by command-line arguments)
DEFAULT_PARALLEL_PROCESSES = 4

def run_anvi_gen_contigs_database(fna_file, colony_output_dir, colony_id):
    """
    Runs anvi-gen-contigs-database.
    """
    logger = logging.getLogger(f'colony_{colony_id}.gen_contigs')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_contigs_db_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Contigs database already created for colony {colony_id}. Skipping.")
        return

    logger.info(f"Creating contigs database for colony {colony_id}...")
    db_file = os.path.join(colony_output_dir, f"{colony_id}_concat.db")
    command = [
        "anvi-gen-contigs-database",
        "-f", fna_file,
        "-o", db_file,
        "-n", colony_id
    ]
    command_str = ' '.join(command)
    logger.info(f"  Running command: {command_str}")
    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.info(f"  Successfully created contigs database: {db_file}")
        logger.debug(f"  Stdout: {result.stdout.decode()}")
        logger.debug(f"  Stderr: {result.stderr.decode()}")
        # Create checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write("Contigs database created")
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error creating contigs database: {db_file}:")
        logger.error(f"  Command: {command_str}")
        logger.error(f"  Return Code: {e.returncode}")
        logger.error(f"  Stdout: {e.stdout.decode()}")
        logger.error(f"  Stderr: {e.stderr.decode()}")

print(run_anvi_gen_contigs_database.__doc__)

def run_ncbi_cogs(colony_output_dir, colony_id):
    """
    Runs anvi-run-ncbi-cogs.
    """
    logger = logging.getLogger(f'colony_{colony_id}.ncbi_cogs')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_ncbi_cogs_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Running NCBI COGs already completed for colony {colony_id}. Skipping.")
        return

    logger.info(f"Running NCBI COGs for colony {colony_id}...")
    db_file = os.path.join(colony_output_dir, f"{colony_id}_concat.db")
    if not os.path.exists(db_file):
        logger.error(f"  Database file not found: {db_file}. Cannot run NCBI COGs.")
        return

    command = [
        "anvi-run-ncbi-cogs",
        "-c", db_file,
        "--num-threads", "4"
    ]
    command_str = ' '.join(command)
    logger.info(f"  Running command: {command_str}")
    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.info(f"  Successfully ran NCBI COGs on {db_file}")
        logger.debug(f"  Stdout: {result.stdout.decode()}")
        logger.debug(f"  Stderr: {result.stderr.decode()}")
        # Create checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write("NCBI COGs completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error running NCBI COGs on {db_file}:")
        logger.error(f"  Command: {command_str}")
        logger.error(f"  Return Code: {e.returncode}")
        logger.error(f"  Stdout: {e.stdout.decode()}")
        logger.error(f"  Stderr: {e.stderr.decode()}")

print(run_ncbi_cogs.__doc__)

def run_hmms(colony_output_dir, colony_id):
    """
    Runs anvi-run-hmms.
    """
    logger = logging.getLogger(f'colony_{colony_id}.hmms')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_hmms_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Running HMMs already completed for colony {colony_id}. Skipping.")
        return

    logger.info(f"Running HMMs for colony {colony_id}...")
    db_file = os.path.join(colony_output_dir, f"{colony_id}_concat.db")
    if not os.path.exists(db_file):
        logger.error(f"  Database file not found: {db_file}. Cannot run HMMs.")
        return

    command = [
        "anvi-run-hmms",
        "-c", db_file,
    ]
    command_str = ' '.join(command)
    logger.info(f"  Running command: {command_str}")
    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.info(f"  Successfully ran HMMs on {db_file}")
        logger.debug(f"  Stdout: {result.stdout.decode()}")
        logger.debug(f"  Stderr: {result.stderr.decode()}")
        # Create checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write("HMMs completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error running HMMs on {db_file}:")
        logger.error(f"  Command: {command_str}")
        logger.error(f"  Return Code: {e.returncode}")
        logger.error(f"  Stdout: {e.stdout.decode()}")
        logger.error(f"  Stderr: {e.stderr.decode()}")

print(run_hmms.__doc__)

def run_scg_taxonomy(colony_output_dir, colony_id):
    """
    Runs anvi-run-scg-taxonomy.
    """
    logger = logging.getLogger(f'colony_{colony_id}.step2')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_step2_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Step 2 (SCG Taxonomy) already completed for colony {colony_id}. Skipping.")
        return

    logger.info(f"Step 2: Running single-copy gene taxonomy for colony {colony_id}...")
    db_file = os.path.join(colony_output_dir, f"{colony_id}_concat.db")
    if not os.path.exists(db_file):
        logger.error(f"  Database file not found: {db_file}. Step 2 cannot be run.")
        return

    command = [
        "anvi-run-scg-taxonomy",
        "-c", db_file,
        "--num-threads", "4"
    ]
    command_str = ' '.join(command)
    logger.info(f"  Running command: {command_str}")
    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.info(f"  Successfully ran taxonomy on {db_file}")
        logger.debug(f"  Stdout: {result.stdout.decode()}")
        logger.debug(f"  Stderr: {result.stderr.decode()}")
        # Create checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write("Step 2 completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error running taxonomy on {db_file}:")
        logger.error(f"  Command: {command_str}")
        logger.error(f"  Return Code: {e.returncode}")
        logger.error(f"  Stdout: {e.stdout.decode()}")
        logger.error(f"  Stderr: {e.stderr.decode()}")

print(run_scg_taxonomy.__doc__)

def generate_profile(colony_output_dir, colony_id, bam_base_dir):
    """
    Generates anvi-profile for each BAM file found.
    """
    logger = logging.getLogger(f'colony_{colony_id}.step3')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_step3_done.txt") # Overall checkpoint removed

    logger.info(f"Step 3: Generating profiles for colony {colony_id}...")
    db_file = os.path.join(colony_output_dir, f"{colony_id}_concat.db")

    bam_search_path = os.path.join(bam_base_dir, "**", colony_id, "**", "*.bam")
    bam_files = glob.glob(bam_search_path, recursive=True)

    if not bam_files:
        logger.warning(f"  No BAM files found for colony {colony_id} in {bam_base_dir}. Skipping profile generation.")
        return

    for bam_file_path in bam_files:
        bam_filename = os.path.splitext(os.path.basename(bam_file_path))[0]
        profile_output_dir = os.path.join(colony_output_dir, f"{bam_filename}-profile")

        command = [
            "anvi-profile",
            "-c", db_file,
            "-o", profile_output_dir,
            "-i", bam_file_path,
            "--profile-SCVs"
        ]
        command_str = ' '.join(command)
        logger.info(f"  Running command: {command_str}")
        try:
            result = subprocess.run(command, check=True, capture_output=True)
            logger.info(f"  Successfully created profile in {profile_output_dir} using {bam_file_path}")
            logger.debug(f"  Stdout: {result.stdout.decode()}")
            logger.debug(f"  Stderr: {result.stderr.decode()}")
        except subprocess.CalledProcessError as e:
            logger.error(f"  Error creating profile for colony {colony_id} using {bam_file_path}:")
            logger.error(f"  Command: {command_str}")
            logger.error(f"  Return Code: {e.returncode}")
            logger.error(f"  Stdout: {e.stdout.decode()}")
            logger.error(f"  Stderr: {e.stderr.decode()}")
            # Optionally, decide if you want to stop processing other BAMs for this colony
            # return

    # Create a general checkpoint for step 3 after processing all BAMs
    with open(checkpoint_file, 'w') as f:
        f.write("Step 3 completed")

print(generate_profile.__doc__)

def merge_profiles_for_colony(colony_output_dir, colony_id):
    """
    Merges all anvi-profiles for a colony.
    """
    logger = logging.getLogger(f'colony_{colony_id}.step4')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_step4_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Step 4 already completed for colony {colony_id}. Skipping.")
        return

    logger.info(f"Step 4: Merging profiles for colony {colony_id}...")
    merged_output_dir = os.path.join(colony_output_dir, f"{colony_id}_merged-profile") # Updated output name
    contigs_db = os.path.join(colony_output_dir, f"{colony_id}_concat.db")
    profile_dirs = glob.glob(os.path.join(colony_output_dir, "*-profile")) # Find all profile directories

    if not os.path.exists(contigs_db):
        logger.warning(f"  Contigs database not found: {contigs_db}. Skipping merge for colony {colony_id}.")
        return
    if not profile_dirs:
        logger.warning(f"  No profile directories found in {colony_output_dir} for colony {colony_id}. Skipping merge.")
        return

    profile_db_paths = [os.path.join(p, "PROFILE.db") for p in profile_dirs]

    command = [
        "anvi-merge",
        *profile_db_paths, # unpack the list of profile DBs
        "-c", contigs_db,
        "-o", merged_output_dir
    ]
    command_str = ' '.join(command)
    logger.info(f"  Running command: {command_str}")
    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.info(f"  Successfully merged profiles for colony {colony_id} in {merged_output_dir}")
        logger.debug(f"  Stdout: {result.stdout.decode()}")
        logger.debug(f"  Stderr: {result.stderr.decode()}")
        # Create checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write("Step 4 completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error merging profiles for colony {colony_id}:")
        logger.error(f"  Command: {command_str}")
        logger.error(f"  Return Code: {e.returncode}")
        logger.error(f"  Stdout: {e.stdout.decode()}")
        logger.error(f"  Stderr: {e.stderr.decode()}")

print(merge_profiles_for_colony.__doc__)

def process_colony(fna_file, colony_id, output_base_dir, bam_base_dir):
    """
    Processes all steps for a single colony.
    """
    colony_output_dir = os.path.join(output_base_dir, f"Colony_{colony_id}")
    os.makedirs(colony_output_dir, exist_ok=True) # Create colony output directory

    log_file = os.path.join(colony_output_dir, f"colony_{colony_id}.log")
    logger = logging.getLogger(f'colony_{colony_id}')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    logger.info(f"Starting processing for colony {colony_id} in directory: {colony_output_dir}")

    run_anvi_gen_contigs_database(fna_file, colony_output_dir, colony_id)
    run_ncbi_cogs(colony_output_dir, colony_id)
    run_hmms(colony_output_dir, colony_id)
    run_scg_taxonomy(colony_output_dir, colony_id)
    generate_profile(colony_output_dir, colony_id, bam_base_dir)
    merge_profiles_for_colony(colony_output_dir, colony_id)

    logger.info(f"Finished processing for colony {colony_id}")
    logging.shutdown()

print(process_colony.__doc__)

def main():
    """
    Main function to orchestrate the Anvi'o workflow.
    """
    parser = argparse.ArgumentParser(description="Anvi'o workflow for processing colonies.")
    parser.add_argument("--fna_dir", required=True, help="Input directory for FNA files.")
    parser.add_argument("--output_dir", required=True, help="Output base directory.")
    parser.add_argument("--bam_base_dir", required=True, help="Base directory for BAM files.")
    parser.add_argument("--processes", type=int, default=DEFAULT_PARALLEL_PROCESSES, help="Number of parallel processes.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    fna_files = glob.glob(os.path.join(args.fna_dir, "*.fna"))
    if not fna_files:
        logging.warning(f"No .fna files found in {args.fna_dir}. Exiting.")
        return

    colony_tasks = []
    for fna_file in fna_files:
        filename_base = os.path.splitext(os.path.basename(fna_file))[0]
        colony_id = filename_base.split("_")[0]
        colony_tasks.append((fna_file, colony_id, args.output_dir, args.bam_base_dir))

    try:
        with multiprocessing.Pool(processes=args.processes) as pool:
            pool.starmap(process_colony, colony_tasks)
    except Exception as e:
        logging.error(f"An error occurred during parallel processing: {e}")
        return

    logging.info("Anvi'o workflow finished for all colonies (parallel).")
    logging.shutdown()

print(main.__doc__)

if __name__ == "__main__":
    main()
