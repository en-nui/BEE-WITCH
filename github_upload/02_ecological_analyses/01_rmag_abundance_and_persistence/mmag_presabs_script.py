import os
import pandas as pd
import csv
import gzip
from collections import defaultdict
import numpy as np
from multiprocessing import Pool

# --- Global File Paths ---
# Directory containing your gzipped per-base BED files
bed_file_dir_global = "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/mosdepth_to_download/"
# Path to your RENAMED_CONTIG_TAXONOMY CSV file
renamed_csv_path_global = "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"
# Base directory where all analysis outputs (per-file CSVs and concatenated CSVs) will be saved
output_base_dir_global = "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/"

# --- Output Headers for Passed Contigs ---
PASSED_HEADERS = [
    'sample_name', 'State', 'Month', 'Caste', 'Colony', 'Rep',
    'MAG_id', 'Genus', 'genome_size', 'type', 'contig_length', 'updated_contig_id',
    'total_contigs_in_MAG_in_sheet', 'observed_contigs_in_bed_file', 'contig_presence_ratio',
    'median_coverage_for_observed_contigs', 'mean_coverage_for_observed_contigs'
]

# --- Output Headers for Filtered Contigs ---
FILTERED_HEADERS = [
    'sample_name', 'State', 'Month', 'Caste', 'Colony', 'Rep',
    'MAG_id', 'Genus', 'genome_size', 'type', 'contig_length', 'updated_contig_id',
    'filter_reason', 'filter_value',
    'total_contigs_in_MAG_in_sheet', 'observed_contigs_in_bed_file',
    'contig_presence_ratio_at_filter', 'median_coverage_at_filter', 'mean_coverage_at_filter'
]

def create_filtered_record(sample_name, state, month, caste, colony, rep,
                           mag_id, genus, genome_size, mag_type, contig_length, updated_contig_id,
                           filter_reason, filter_value,
                           total_contigs_in_MAG_in_sheet, observed_contigs_in_bed_file,
                           contig_presence_ratio_at_filter, median_coverage_at_filter, mean_coverage_at_filter):
    """Helper to create a dictionary record for filtered contigs."""
    return {
        'sample_name': sample_name, 'State': state, 'Month': month, 'Caste': caste, 'Colony': colony, 'Rep': rep,
        'MAG_id': mag_id, 'Genus': genus, 'genome_size': genome_size, 'type': mag_type,
        'contig_length': contig_length, 'updated_contig_id': updated_contig_id,
        'filter_reason': filter_reason, 'filter_value': filter_value,
        'total_contigs_in_MAG_in_sheet': total_contigs_in_MAG_in_sheet,
        'observed_contigs_in_bed_file': observed_contigs_in_bed_file,
        'contig_presence_ratio_at_filter': contig_presence_ratio_at_filter,
        'median_coverage_at_filter': median_coverage_at_filter,
        'mean_coverage_at_filter': mean_coverage_at_filter
    }

def process_single_bed_file(
    bed_file_name: str,
    bed_file_dir: str,
    renamed_data_maps: dict,
    output_base_dir: str
) -> tuple:
    """
    Processes a single gzipped BED file, applies MAG filters, and returns
    lists of passed and filtered contig records.
    This function is designed to be called by a multiprocessing Pool.
    """

    passed_contig_records_for_sample = []
    filtered_contig_records_for_sample = []

    current_bed_path = os.path.join(bed_file_dir, bed_file_name)

    # Unpack renamed_data_maps
    contig_to_metadata_map = renamed_data_maps['contig_to_metadata_map']
    mag_to_all_contigs_in_sheet = renamed_data_maps['mag_to_all_contigs_in_sheet']

    # --- Extract Sample Metadata from filename ---
    state, month, caste, colony, rep, sample_name = None, None, None, None, None, None
    try:
        base_name = bed_file_name.replace(".per-base.bed.gz", "")
        parts = base_name.split('_')
        if len(parts) >= 5:
            state = parts[0]
            month = parts[1]
            caste = parts[2]
            colony = int(parts[3])
            rep = int(parts[4])
            sample_name = "_".join(parts[:5])
        else:
            print(f"[Worker] Warning: Filename '{bed_file_name}' does not match expected format. Skipping metadata extraction.")
    except ValueError:
        print(f"[Worker] Warning: Could not parse Colony or Rep as integer from '{bed_file_name}'. Skipping metadata extraction.")
    except Exception as e:
        print(f"[Worker] Error extracting metadata from filename '{bed_file_name}': {e}. Skipping metadata extraction.")

    print(f"[Worker] Processing: {bed_file_name}")

    mag_coverage_data = defaultdict(lambda: {'coverage_segments': []}) # Stores (coverage, length) for each segment
    mag_observed_contigs = defaultdict(set) # Tracks contigs observed in this BED file, by MAG_id

    try:
        with gzip.open(current_bed_path, 'rt') as bed_in:
            for line_num, line in enumerate(bed_in, 1):
                parts = line.strip().split('\t')
                if len(parts) != 4:
                    continue

                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    coverage = float(parts[3])
                except ValueError:
                    continue

                contig_length_in_bed_segment = end - start
                if contig_length_in_bed_segment <= 0:
                    continue

                # --- Map chrom to RENAMED and handle immediate filters ---
                if chrom not in contig_to_metadata_map:
                    # Contig not found in RENAMED: immediate filter for this specific contig
                    # Added to filtered list as it's an mmag if it were in the sheet.
                    # This case handles bed file contigs not present in RENAMED sheet at all.
                    filtered_contig_records_for_sample.append(create_filtered_record(
                        sample_name, state, month, caste, colony, rep,
                        'N/A', 'N/A', 'N/A', 'N/A', 'N/A', chrom,
                        "contig_not_in_renamed", 'N/A', 0, 0, 0.0, 0.0, 0.0
                    ))
                    print(f"[Worker] Warning: Contig '{chrom}' from {bed_file_name} not found in RENAMED_CONTIG_TAXONOMY. Ignoring this contig.")
                    continue

                contig_metadata = contig_to_metadata_map[chrom]
                contig_type = contig_metadata.get('type')

                # NEW LOGIC: Only process contigs where type == 'mmag'. Skip all others.
                if contig_type != 'mmag':
                    # Do not add to filtered_contig_records_for_sample, simply skip.
                    continue

                # This is an mmag contig, proceed with its MAG data accumulation
                mag_id = contig_metadata.get('MAG_id')
                if mag_id is None:
                    # Contig has no MAG_id in RENAMED: immediate filter
                    filtered_contig_records_for_sample.append(create_filtered_record(
                        sample_name, state, month, caste, colony, rep,
                        'N/A', contig_metadata.get('Genus', 'N/A'), contig_metadata.get('genome_size', 'N/A'),
                        contig_type, contig_metadata.get('contig_length', 'N/A'), chrom,
                        "contig_has_no_MAG_id", 'N/A', 0, 0, 0.0, 0.0, 0.0
                    ))
                    print(f"[Worker] Warning: Contig '{chrom}' has no 'MAG_id' in RENAMED. Ignoring this contig.")
                    continue

                mag_observed_contigs[mag_id].add(chrom)
                mag_coverage_data[mag_id]['coverage_segments'].append((coverage, contig_length_in_bed_segment))

    except FileNotFoundError:
        print(f"[Worker] Error: BED file not found at {current_bed_path}. Skipping.")
        return [], []
    except Exception as e:
        print(f"[Worker] An unexpected error occurred processing {bed_file_name}: {e}. Skipping.")
        return [], []

    # --- Apply MAG-level Filters for mmag-types ---
    passed_mags_data = {}

    for mag_id, observed_contigs_set in mag_observed_contigs.items():
        representative_contig_id = next(iter(observed_contigs_set), None)
        if representative_contig_id is None: continue

        mag_level_metadata = contig_to_metadata_map.get(representative_contig_id, {})
        mag_genome_size = mag_level_metadata.get('genome_size')
        mag_genus = mag_level_metadata.get('Genus', 'N/A')
        mag_type = mag_level_metadata.get('type') # This should always be 'mmag' here now

        current_total_mag_contigs_in_sheet = len(mag_to_all_contigs_in_sheet.get(mag_id, set()))
        current_observed_contigs_in_bed_file = len(observed_contigs_set)

        current_contig_presence_ratio = 0.0
        if current_total_mag_contigs_in_sheet > 0:
            current_contig_presence_ratio = current_observed_contigs_in_bed_file / current_total_mag_contigs_in_sheet

        current_median_coverage = 0.0
        current_mean_coverage = 0.0

        current_mag_cov_segments = mag_coverage_data[mag_id]['coverage_segments']

        # Calculate mean and median for the MAG
        total_bases_covered = sum(seg[1] for seg in current_mag_cov_segments)
        sum_coverage_values = sum(seg[0] * seg[1] for seg in current_mag_cov_segments)

        if total_bases_covered > 0:
            current_mean_coverage = sum_coverage_values / total_bases_covered

            coverage_values_for_median = [seg[0] for seg in current_mag_cov_segments]
            segment_lengths_for_median = [seg[1] for seg in current_mag_cov_segments]

            try:
                expanded_coverages = np.repeat(coverage_values_for_median, segment_lengths_for_median)
                current_median_coverage = np.median(expanded_coverages)
            except MemoryError:
                print(f"[Worker] Warning: MemoryError calculating median for MAG_id {mag_id} in {bed_file_name}. Treating as failed median threshold.")
                current_median_coverage = -1 # Forces filter
            except Exception as e:
                print(f"[Worker] Error calculating median for MAG_id {mag_id} in {bed_file_name}: {e}. Setting median to -1.")
                current_median_coverage = -1

        filter_reason = None
        filter_value = None

        # Filter 1: Genome Size
        if mag_genome_size is not None and isinstance(mag_genome_size, (int, float)) and mag_genome_size > 10_000_000:
            filter_reason = "MAG_genome_size_too_large"
            filter_value = mag_genome_size
        # Filter 2: Contig Presence Ratio
        elif current_contig_presence_ratio < 0.25:
            filter_reason = "MAG_contig_presence_ratio_low"
            filter_value = current_contig_presence_ratio
        # Filter 3: Coverage Thresholds (NEW LOGIC)
        # Keep if: median_coverage >= 5
        # OR Keep if: (median_coverage < 5 AND median_coverage >= 1) AND (mean_coverage > 1)
        # Otherwise, remove.
        elif not (current_median_coverage >= 1 or \
                  (current_median_coverage < 5 and current_median_coverage >= 1 and current_mean_coverage > 1)):
            filter_reason = "MAG_coverage_threshold_not_met"
            filter_value = f"Median:{current_median_coverage:.2f}, Mean:{current_mean_coverage:.2f}"

        if filter_reason:
            # MAG was filtered out: record all its observed contigs for the filtered report
            for chrom in observed_contigs_set:
                contig_metadata_for_filtered = contig_to_metadata_map.get(chrom, {})
                filtered_contig_records_for_sample.append(create_filtered_record(
                    sample_name, state, month, caste, colony, rep,
                    mag_id, contig_metadata_for_filtered.get('Genus', 'N/A'),
                    contig_metadata_for_filtered.get('genome_size', 'N/A'),
                    contig_metadata_for_filtered.get('type'),
                    contig_metadata_for_filtered.get('contig_length'), chrom,
                    filter_reason, filter_value,
                    current_total_mag_contigs_in_sheet, current_observed_contigs_in_bed_file,
                    current_contig_presence_ratio, current_median_coverage, current_mean_coverage
                ))
        else:
            # MAG passed all filters, store its calculated metrics
            passed_mags_data[mag_id] = {
                'Genus': mag_genus, 'genome_size': mag_genome_size, 'type': mag_type,
                'total_contigs_in_MAG_in_sheet': current_total_mag_contigs_in_sheet,
                'observed_contigs_in_bed_file': current_observed_contigs_in_bed_file,
                'contig_presence_ratio': current_contig_presence_ratio,
                'median_coverage_for_observed_contigs': current_median_coverage,
                'mean_coverage_for_observed_contigs': current_mean_coverage
            }

    # --- Populate passed_contig_records_for_sample ---
    for mag_id, mag_metrics in passed_mags_data.items():
        for chrom in mag_observed_contigs[mag_id]: # Iterate over only observed contigs that passed
            contig_metadata = contig_to_metadata_map.get(chrom)
            if contig_metadata is None: # Should not happen given prior checks, but for safety
                continue

            record = {
                'sample_name': sample_name, 'State': state, 'Month': month, 'Caste': caste, 'Colony': colony, 'Rep': rep,
                'MAG_id': mag_id, 'Genus': mag_metrics['Genus'], 'genome_size': mag_metrics['genome_size'], 'type': mag_metrics['type'],
                'contig_length': contig_metadata.get('contig_length'), 'updated_contig_id': chrom,
                'total_contigs_in_MAG_in_sheet': mag_metrics['total_contigs_in_MAG_in_sheet'],
                'observed_contigs_in_bed_file': mag_metrics['observed_contigs_in_bed_file'],
                'contig_presence_ratio': mag_metrics['contig_presence_ratio'],
                'median_coverage_for_observed_contigs': mag_metrics['median_coverage_for_observed_contigs'],
                'mean_coverage_for_observed_contigs': mag_metrics['mean_coverage_for_observed_contigs']
            }
            passed_contig_records_for_sample.append(record)

    # --- Write per-file output CSVs ---
    if passed_contig_records_for_sample:
        per_file_output_name = f"{sample_name}.filtered_mmag_contigs.csv" # Renamed for clarity
        per_file_output_path = os.path.join(output_base_dir, per_file_output_name)

        per_file_df = pd.DataFrame(passed_contig_records_for_sample, columns=PASSED_HEADERS)
        per_file_df.to_csv(per_file_output_path, index=False)
        print(f"[Worker] Generated per-file PASSED CSV: {per_file_output_path} ({len(passed_contig_records_for_sample)} contig records)")
    else:
        print(f"[Worker] No mmag contigs passed filters for {bed_file_name}. No per-file PASSED CSV generated.")

    if filtered_contig_records_for_sample:
        per_file_filtered_output_name = f"{sample_name}.removed_mmag_contigs.csv" # Renamed for clarity
        per_file_filtered_output_path = os.path.join(output_base_dir, per_file_filtered_output_name)
        per_file_filtered_df = pd.DataFrame(filtered_contig_records_for_sample, columns=FILTERED_HEADERS)
        per_file_filtered_df.to_csv(per_file_filtered_output_path, index=False)
        print(f"[Worker] Generated per-file REMOVED CSV: {per_file_filtered_output_path} ({len(filtered_contig_records_for_sample)} contig records)")
    else:
        print(f"[Worker] No mmag contigs were removed by filters for {bed_file_name}. No per-file REMOVED CSV generated.")

    # Return the records to the main process for concatenation
    return passed_contig_records_for_sample, filtered_contig_records_for_sample

def main_analysis_parallel(
    bed_file_dir: str,
    renamed_csv_path: str,
    output_base_dir: str,
    num_cpus: int = 4
):
    """
    Main function to orchestrate parallel processing of BED files.
    """
    print(f"Main Process: Starting parallel analysis with {num_cpus} CPUs...")
    os.makedirs(output_base_dir, exist_ok=True)

    # --- Load RENAMED_CONTIG_TAXONOMY.csv (done once in main process) ---
    try:
        renamed_df = pd.read_csv(renamed_csv_path)
        contig_to_metadata_map = renamed_df.set_index('updated_contig_id').to_dict('index')
        mag_to_all_contigs_in_sheet = defaultdict(set) # Correct variable name
        for _, row in renamed_df.iterrows():
            mag_to_all_contigs_in_sheet[row['MAG_id']].add(row['updated_contig_id'])
        print(f"Main Process: Loaded RENAMED_CONTIG_TAXONOMY.csv with {len(renamed_df)} entries.")

    except FileNotFoundError:
        print(f"Main Process Error: RENAMED_CONTIG_TAXONOMY.csv not found at {renamed_csv_path}")
        return
    except KeyError as e:
        print(f"Main Process Error: Missing expected column in RENAMED_CONTIG_TAXONOMY.csv: {e}")
        print("Please ensure 'updated_contig_id', 'MAG_id', 'type', 'genome_size', 'contig_length', 'Genus' columns exist.")
        return
    except Exception as e:
        print(f"Main Process: An unexpected error occurred loading RENAMED CSV: {e}")
        return

    bed_files = [f for f in os.listdir(bed_file_dir) if f.endswith('.per-base.bed.gz')]
    if not bed_files:
        print(f"Main Process: No *.per-base.bed.gz files found in {bed_file_dir}. Exiting.")
        return

    print(f"Main Process: Found {len(bed_files)} BED files to process.")

    # Package necessary data into a dict to pass to worker processes
    # These will be copied to each worker process
    renamed_data_maps = {
        'contig_to_metadata_map': contig_to_metadata_map,
        'mag_to_all_contigs_in_sheet': mag_to_all_contigs_in_sheet # Corrected variable name here
    }

    # Prepare arguments for multiprocessing map
    # Each item in `tasks` is a tuple of individual arguments
    tasks = [
        (bf, bed_file_dir, renamed_data_maps, output_base_dir)
        for bf in bed_files
    ]

    all_concatenated_passed_results = []
    all_concatenated_filtered_results = []

    # Use multiprocessing Pool to parallelize the processing
    with Pool(processes=num_cpus) as pool:
        # starmap() unpacks each tuple from 'tasks' as separate arguments to 'process_single_bed_file'
        results = pool.starmap(process_single_bed_file, tasks)

    # Collect results from all processes
    for passed_recs, filtered_recs in results:
        all_concatenated_passed_results.extend(passed_recs)
        all_concatenated_filtered_results.extend(filtered_recs)

    # --- Write concatenated passed results ---
    if all_concatenated_passed_results:
        concatenated_passed_output_name = "all_filtered_mmag_contigs_concatenated.csv"
        concatenated_passed_output_path = os.path.join(output_base_dir, concatenated_passed_output_name)

        concatenated_passed_df = pd.DataFrame(all_concatenated_passed_results, columns=PASSED_HEADERS)
        concatenated_passed_df.to_csv(concatenated_passed_output_path, index=False)
        print(f"\nMain Process: All passed results saved to: {concatenated_passed_output_path} ({len(all_concatenated_passed_results)} total contig records)")
    else:
        print("\nMain Process: No mmag contigs passed filters across all BED files. No concatenated PASSED CSV generated.")

    # --- Write concatenated filtered results ---
    if all_concatenated_filtered_results:
        concatenated_filtered_output_name = "all_removed_mmag_contigs_concatenated.csv"
        concatenated_filtered_output_path = os.path.join(output_base_dir, concatenated_filtered_output_name)

        concatenated_filtered_df = pd.DataFrame(all_concatenated_filtered_results, columns=FILTERED_HEADERS)
        concatenated_filtered_df.to_csv(concatenated_filtered_output_path, index=False)
        print(f"\nMain Process: All removed mmag contigs saved to: {concatenated_filtered_output_path} ({len(all_concatenated_filtered_results)} total contig records)")
    else:
        print("\nMain Process: No mmag contigs were removed across all BED files. No concatenated REMOVED CSV generated.")

    print("\nMain Process: Analysis completed.")

if __name__ == "__main__":
    # Ensure the output directory exists before starting parallel processes
    os.makedirs(output_base_dir_global, exist_ok=True)

    # Call the main parallel analysis function
    main_analysis_parallel(
        bed_file_dir=bed_file_dir_global,
        renamed_csv_path=renamed_csv_path_global,
        output_base_dir=output_base_dir_global,
        num_cpus=4 # Set the number of CPUs here
    )
