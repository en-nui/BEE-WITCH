import os
import subprocess
import logging
#i ran this in my mamba env with binette. 
logging.basicConfig(
    filename="binette_pipeline.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
base_metabat_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/07a_mMAG_binning"
base_vamb_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/07b_mMAG_binning_vamb"
base_assembly_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/02_Assemblies"
output_base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/10_binette_results"

def find_and_process_bin_dirs(base_metabat_dir, base_vamb_dir, base_assembly_dir, output_base_dir):
    """
    Finds metbat bin and corresponding vamb/assembly paths. Processes with binette.
    
    Args:
        base_metabat_dir (str): Base directory for metabat bins.
        base_vamb_dir (str): Base directory for vamb bins.
        base_assembly_dir (str): Base directory for assemblies.
        output_base_dir (str): Base directory for binette output.
    """
for root, dirs, files in os.walk(base_metabat_dir):
    metabat_bins = [f for f in files if f.endswith(".fa")]
    if not metabat_bins:
        continue
    #this checks to see if the current directory contains *.fa files in metabat

    relpath = os.path.relpath(root, base_metabat_dir)
    
    #find the paths
    vamb_dir = os.path.join(base_vamb_dir, relpath, "cluster_fasta")
    assembly_path = os.path.join(base_assembly_dir, relpath, "final.contigs.fa")
    output_dir = os.path.join(output_base_dir, relpath)
    
    #logging paths!!!!!!!
    logging.info(f"Processing paths for sample: {relpath}")
    logging.info(f"Metabat dir: {root}")
    logging.info(f"Vamb dir: {vamb_dir}")
    logging.info(f"Assembly path: {assembly_path}")
    logging.info(f"Output dir: {output_dir}")

    #find corresponding vamb directory and assembly file
    if not os.path.exists(vamb_dir):
        logging.warning(f"Vamb directory does not exist for {relpath}. Skipping.")
        continue

    vamb_bins = [f for f in os.listdir(vamb_dir) if f.endswith(".fna")]
    if not vamb_bins:
        logging.warning(f"No vamb bins found in {vamb_dir}. Skipping.")
        continue
    if not os.path.exists(assembly_path):
        logging.warning(f"Assembly file does not exist for {relpath}. Skipping.")
        continue

    #check to see if output directory has already been processed
    checkpoint_file = os.path.join(output_dir, "binette_complete.txt")
    if os.path.exists(checkpoint_file):
        logging.info(f"Checkpoint file found. Skipping already processed sample: {relpath}")
        continue

    #create output directory
    try:
        command = [
            "binette",
            "--bin_dirs", root, vamb_dir,
            "--contigs", assembly_path,
            "--outdir", output_dir
        ]

        subprocess.run(command, check=True)
        logging.info(f"Successfully processed with binette: {relpath}")

    #create checkpoint
    
        with open(checkpoint_file, "w") as f:
            f.write("Binette processing complete.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running binette for {relpath}: {e}")
    except Exception as e:
        logging.error(f"Unexpected error for {relpath}: {e}")

if __name__ == "__main__":
    find_and_process_bin_dirs(
        base_metabat_dir,
        base_vamb_dir,
        base_assembly_dir,
        output_base_dir,
    )
