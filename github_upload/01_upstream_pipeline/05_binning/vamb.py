import os
import subprocess
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)

# Define directory paths
mapping_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/03_Mapping"
assembly_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/02_Assemblies"
output_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/07b_mMAG_binning_vamb"

# Function to find deepest directories
def get_deepest_dirs(base_dir):
    return [
        str(path) for path in Path(base_dir).rglob("*")
        if path.is_dir() and not any(child.is_dir() for child in path.iterdir())
    ]

# Function to process directories and check for BAM/BAI files
def process_directory(subdir):
    logging.info(f"Checking directory: {subdir}")
    bam_files = [f for f in Path(subdir).glob("*.bam")]
    bai_files = [f for f in Path(subdir).glob("*.bam.bai")]

    bam_bai_pairs = [
        bam for bam in bam_files if Path(f"{bam}.bai").exists()
    ]
    has_files = bool(bam_bai_pairs)
    if not has_files:
        logging.warning(f"Missing BAM/BAI pair in {subdir}")
    return has_files

# Function to group subdirectories by parent
def group_by_parent(subdirs):
    grouped = {}
    for subdir in subdirs:
        parent = str(Path(subdir).parent)
        grouped.setdefault(parent, []).append(subdir)
    return grouped

# Function to check all subdirectories under each parent
def validate_parent_dirs(grouped_dirs):
    valid_parents = []
    for parent, subdirs in grouped_dirs.items():
        logging.info(f"Validating parent directory: {parent}")
        if all(process_directory(subdir) for subdir in subdirs):
            valid_parents.append(parent)
        else:
            logging.warning(f"Skipping parent directory {parent}: Incomplete BAM/BAI pairs")
    return valid_parents

# Function to find coassembly file
def find_coassembly(parent_dir):
    hierarchy = "/".join(Path(parent_dir).parts[-4:])  # Adjust based on directory depth
    return os.path.join(assembly_dir, hierarchy, "final.contigs.fa")

# Function to check if VAMB has already been run
def is_vamb_complete(output_path):
    expected_file = os.path.join(output_path, "vamb_output.1.fa")
    return os.path.exists(expected_file)

# Function to run VAMB
def run_vamb(mapping_subdirs, coassembly, output_path):
    if is_vamb_complete(output_path):
        logging.info(f"Skipping VAMB: Output already exists at {output_path}")
        return

    bam_files = [str(bam) for subdir in mapping_subdirs for bam in Path(subdir).glob("*.bam")]
    if not bam_files:
        logging.warning(f"No BAM files found for coassembly {coassembly}")
        return

    # Do not create the directory manually
    cmd = ["vamb", "bin", "default", "--outdir", output_path, "--fasta", coassembly, "--bamfiles", *bam_files]
    logging.info(f"Running VAMB: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"VAMB failed for {coassembly}: {e}")

# Main execution logic
if __name__ == "__main__":
    logging.info("Scanning for valid directories...")
    all_dirs = get_deepest_dirs(mapping_dir)
    logging.info(f"Found directories: {all_dirs}")

    # Group directories by parent
    grouped_dirs = group_by_parent(all_dirs)
    logging.info(f"Grouped directories by parent: {grouped_dirs}")

    # Validate parent directories
    valid_parents = validate_parent_dirs(grouped_dirs)

    # Logging results
    logging.info(f"Valid parent directories: {valid_parents}")
    
    # Ensure the base output directory exists (but not the final subdirectories)
    os.makedirs(output_dir, exist_ok=True)

    # Create necessary parent directories but not the final directory
    for parent in valid_parents:
        # Create the directory structure for the parent, but not the deepest subdir (e.g., /100)
        relative_path = Path(parent).relative_to(mapping_dir)  # Relative path under the mapping directory
        partial_output_path = os.path.join(output_dir, *relative_path.parts[:-1])  # Exclude the final subdir

        os.makedirs(partial_output_path, exist_ok=True)
    
    # Write log file at the base directory
    with open(os.path.join(output_dir, "bam_check_log.txt"), "w") as log:
        for parent in valid_parents:
            log.write(f"{parent}: All BAM and BAI files found\n")
        for parent, subdirs in grouped_dirs.items():
            if parent not in valid_parents:
                log.write(f"{parent}: Missing BAM or BAI files in some subdirectories\n")

    # Process each valid parent directory
    logging.info("Starting VAMB binning...")
    with ProcessPoolExecutor(max_workers=24) as executor:
        futures = []
        for parent in valid_parents:
            coassembly_path = find_coassembly(parent)
            print(coassembly_path)
            if os.path.exists(coassembly_path):
                # Generate the corresponding output path, matching the assembly path structure
                relative_path = Path(parent).relative_to(mapping_dir)  # Relative path under the mapping directory
                output_path = os.path.join(output_dir, *relative_path.parts)
                print(output_path)
                futures.append(executor.submit(run_vamb, grouped_dirs[parent], coassembly_path, output_path))
            else:
                logging.warning(f"Coassembly file not found for {parent}")

        for future in futures:
            future.result()

    logging.info("VAMB binning complete.")

