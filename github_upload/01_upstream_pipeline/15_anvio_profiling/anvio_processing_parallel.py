def generate_profile(colony_output_dir, colony_id):
    """
    Generates anvi-profile.
    """
    logger = logging.getLogger('step3')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_step3_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Step 3 already completed for colony {colony_id}. Skipping.")
        return

    logger.info(f"Step 3: Generating profile for colony {colony_id}...")
    db_file = os.path.join(colony_output_dir, f"{colony_id}_concat.db") # Assuming DB named like this
    profile_output_dir = os.path.join(colony_output_dir, f"{colony_id}_concat-profile")
    bam_file_path = os.path.join(
        BAM_BASE_DIR, "IL", "June", "Worker", colony_id, "1", f"IL_June_Worker_{colony_id}_1.bam"
    )

    if not os.path.exists(bam_file_path):
        logger.warning(f"  BAM file not found: {bam_file_path}. Skipping profile generation for colony {colony_id}.")
        return

    command = [
        "anvi-profile",
        "-c", db_file,
        "-o", profile_output_dir,
        "-i", bam_file_path
    ]
    command_str = ' '.join(command)
    logger.info(f"  Running command: {command_str}")
    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.info(f"  Successfully created profile in {profile_output_dir}")
        logger.debug(f"  Stdout: {result.stdout.decode()}")
        logger.debug(f"  Stderr: {result.stderr.decode()}")
        # Create checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write("Step 3 completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error creating profile for colony {colony_id} using {bam_file_path}:")
        logger.error(f"  Command: {command_str}")
        logger.error(f"  Return Code: {e.returncode}")
        logger.error(f"  Stdout: {e.stdout.decode()}")
        logger.error(f"  Stderr: {e.stderr.decode()}")

def merge_profiles_for_colony(colony_output_dir, colony_id):
    """
    Merges anvi-profiles for a colony.
    """
    logger = logging.getLogger('step4')
    checkpoint_file = os.path.join(colony_output_dir, f"{colony_id}_step4_done.txt")
    if os.path.exists(checkpoint_file):
        logger.info(f"Step 4 already completed for colony {colony_id}. Skipping.")
        return

    logger.info(f"Step 4: Merging profiles for colony {colony_id}...")
    merged_output_dir = os.path.join(colony_output_dir, f"{colony_id}_concat_merged-profiles")
    contigs_db = os.path.join(colony_output_dir, f"{colony_id}_concat.db")
    profile_dirs = glob.glob(os.path.join(colony_output_dir, f"{colony_id}_concat-profile")) # Expecting only one

    if not os.path.exists(contigs_db):
        logger.warning(f"  Contigs database not found: {contigs_db}. Skipping merge for colony {colony_id}.")
        return
    if not profile_dirs:
        logger.warning(f"  No profile directories found in {colony_output_dir} for colony {colony_id}. Skipping merge.")
        return

    profile_db_paths = [os.path.join(profile_dirs[0], "PROFILE.db")] # Assuming only one profile dir
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


def process_colony(fna_file, colony_id, output_base_dir):
    """
    Processes all steps for a single colony.
    """
    colony_output_dir = os.path.join(output_base_dir, f"Colony_{colony_id}")
    os.makedirs(colony_output_dir, exist_ok=True) # Create colony output directory

    log_file = os.path.join(colony_output_dir, f"colony_{colony_id}.log")
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) # Set root logger level to DEBUG to capture all levels
    fh = logging.FileHandler(log_file, mode='w') # Log to file in colony dir
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s') # Include logger name
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    console_handler = logging.StreamHandler() # Also log to console
    console_handler.setLevel(logging.INFO) # But console only for INFO and above
    console_formatter = logging.Formatter('%(levelname)s - %(message)s') # Simpler format for console
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)


    logger.info(f"Starting processing for colony {colony_id} in directory: {colony_output_dir}")

    run_anvi_gen_contigs_database(fna_file, colony_output_dir, colony_id)
    run_scg_taxonomy(colony_output_dir, colony_id)
    generate_profile(colony_output_dir, colony_id)
    merge_profiles_for_colony(colony_output_dir, colony_id)

    logger.info(f"Finished processing for colony {colony_id}")
    logging.shutdown() # Close handlers for next colony processing if any


def main():
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True) # Ensure base output directory exists

    fna_files = glob.glob(os.path.join(FNA_INPUT_DIR, "*.fna"))
    if not fna_files:
        logging.warning(f"No .fna files found in {FNA_INPUT_DIR}. Exiting.")
        return

    for fna_file in fna_files:
        filename_base = os.path.splitext(os.path.basename(fna_file))[0]
        colony_id = filename_base.split("_")[0] # e.g., "100" from "100_concat.fna"
        process_colony(fna_file, colony_id, OUTPUT_BASE_DIR)


    logging.info("Anvi'o workflow finished for all colonies.")

if __name__ == "__main__":
    main()

