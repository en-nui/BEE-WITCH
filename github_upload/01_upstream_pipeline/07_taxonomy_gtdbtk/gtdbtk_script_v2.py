#functionally similar to first script except this includes some error handling and logging 
#also, instead of running gtdbtk on everything within a bini, it new checks to make sure that the relative 'completeness' and 'contamination' scores are above defined thresholds before running the function

import os
import logging
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import shutil
def setup_logging(log_file):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def filter_bins(tsv_path, min_completeness=50, max_contamination=10):
    #parses tsv file within base_path to filter bins based on completeness and contamination scores as generated via checkM@

    logging.info(f"Parsing and filtering bins from {tsv_path}...")
    df = pd.read_csv(tsv_path, sep="\t")

    filtered_bins = df[(df['completeness'] >= min_completeness) & (df['contamination'] <= max_contamination)]
    return filtered_bins['bin_id'].tolist(), len(df), len(filtered_bins)

def prepare_input_bins(bin_ids, final_bins_dir, temp_input_dir):
    #copy filtered bins from filter_bins from the final_bins_dir to the temp_dir

    logging.info(f"Preparing input bins in {temp_input_dir}...")
    os.makedirs(temp_input_dir, exist_ok=True)
    for bin_id in bin_ids:
        bin_file = os.path.join(final_bins_dir, f"bin_{bin_id}.fa")
        if os.path.exists(bin_file):
            destination = os.path.join(temp_input_dir, f"bin_{bin_id}.fna")
            shutil.copy(bin_file, destination)
            logging.info(f"Copied {bin_file} to {destination}")
        else:
            logging.warning(f"Bin file not found: {bin_file}")
            



def run_gtdbtk(bin_directory, output_directory, cpus):
    """
    Runs the GTDB-Tk workflow for a given bin directory.
    """
    try:
        logging.info(f"Starting GTDB-Tk for: {bin_directory}")
        os.makedirs(output_directory, exist_ok=True)  # Ensure the output directory exists
        command = [
            "gtdbtk", "classify_wf",
            "--genome_dir", bin_directory,
            "--out_dir", output_directory,
            "--cpus", str(cpus)
        ]

        subprocess.run(command, check=True)
        logging.info(f"Successfully completed: {bin_directory}")
        return bin_directory, "success"
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to complete: {bin_directory}. Error: {e}")
        return bin_directory, "failure"


def main():
    base_path = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/10_binette_results"
    output_base = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/11_gtdbtk_output"
    log_file = os.path.join(output_base, "gtdbtk_run.log")
    cpus = 4
    max_tasks = 4

    os.makedirs(output_base, exist_ok=True)

    setup_logging(log_file)

    #iterate through directoryes under base_path
    for root, dirs, files in os.walk(base_path):
        if "final_bins_quality_reports.tsv" in files:
            sample_path = root
            logging.info(f"Processing sample: {sample_path}")

            tsv_path = os.path.join(sample_path, "final_bins_quality_reports.tsv")
            final_bins_dir = os.path.join(sample_path, "final_bins")
            temp_input_dir = os.path.join(output_base, os.path.relpath(sample_path, base_path), "gtdbtk_input")
            output_directory = os.path.join(output_base, os.path.relpath(sample_path, base_path), "gtdbtk_output")
            summary_file = os.path.join(output_base, "summary.log")
    
            bin_ids, total_bins, filtered_bins_count = filter_bins(tsv_path)
            prepare_input_bins(bin_ids, final_bins_dir, temp_input_dir)

            results = []
            try:
                with ProcessPoolExecutor(max_tasks) as executor:
                    future = executor.submit(run_gtdbtk, temp_input_dir, output_directory, cpus)
                    results.append(future.result())
            finally:
                logging.info(f"Cleaning up temp input directory: {temp_input_dir}")
                shutil.rmtree(temp_input_dir, ignore_errors=True)

    #logging
    failed_bins = total_bins - filtered_bins_count
    with open(summary_file, "a") as sf:
        sf.write(f"Sample: {sample_path}\n")
        sf.write(f"Total bins: {total_bins}\n")
        sf.write(f"Bins passing QC: {filtered_bins_count}\n")
        sf.write(f"Bins ignored: {failed_bins}\n")
        sf.write(f"Results: {results}\n\n")


if __name__ == "__main__":
    main()


