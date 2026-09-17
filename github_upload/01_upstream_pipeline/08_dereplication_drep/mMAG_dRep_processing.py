import os
import logging
import pandas as pd 
import shutil
from collections import defaultdict

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
    logging.info(f"Parsing and filtering bins from {tsv_path}")
    df = pd.read_csv(tsv_path, sep='\t')

    filtered_bins = df[(df['completeness'] >= min_completeness) & (df['contamination'] <= max_contamination)]

    return filtered_bins['bin_id'].tolist(), len(df), len(filtered_bins)



def format_bin_filename(bin_id, source_path):
    rel_path = os.path.relpath(source_path, start="/N/project/NewtonLab/robinch/projects/BEE-WITCH/10_binette_results")

    formatted_path = rel_path.replace(os.sep,"_")
    return f"{formatted_path}_bin_{bin_id}.fna"

def copy_filtered_bins(bin_ids, source_dir, destination_dir, bin_log):
    os.makedirs(destination_dir, exist_ok=True)

    for bin_id in bin_ids:
        bin_file = os.path.join(source_dir, f"bin_{bin_id}.fa")
        if os.path.exists(bin_file):
            filename = format_bin_filename(bin_id, source_dir)
            destination = os.path.join(destination_dir, filename)
            shutil.copy(bin_file, destination)
            logging.info(f"Copied {bin_file} to {destination}")
            bin_log.append((bin_file, destination))
        else:
            logging.warning(f"Bin file not found: {bin_file}")

def main():
    base_path = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/10_binette_results"
    output_base = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/14_mMAG_dRep"
    temp_directory = os.path.join(output_base, "tmp_directory")
    log_file = os.path.join(output_base, "mMAG_dRep.log")
    summary_file = os.path.join(output_base, "summary.log")

    os.makedirs(output_base, exist_ok=True)
    setup_logging(log_file)

    total_bins, filtered_bins_count = 0,0
    bin_log = []

    for root, dirs, files in os.walk(base_path):
        if "final_bins_quality_reports.tsv" in files:
            tsv_path = os.path.join(root, "final_bins_quality_reports.tsv")
            final_bins_dir = os.path.join(root, "final_bins")
        
            bin_ids, sample_total, sample_filtered = filter_bins(tsv_path)
            copy_filtered_bins(bin_ids, final_bins_dir, temp_directory, bin_log)

            total_bins += sample_total
            filtered_bins_count += sample_filtered

    with open(summary_file, "a") as sf:
        sf.write(f"Total bins: {total_bins}\n")
        sf.write(f"Bins passing QC: {filtered_bins_count}\n")
        sf.write("Copied Files:\n")
        for src, dest in bin_log:
            sf.write(f"{src} -> {dest}\n")

    logging.info("Processing completed.")

if __name__ == "__main__":
    main()
            

