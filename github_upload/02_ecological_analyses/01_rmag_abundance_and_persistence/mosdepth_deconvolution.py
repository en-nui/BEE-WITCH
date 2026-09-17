import os
import pandas as pd
import csv

def analyze_mosdepth_summaries(
    summary_files_dir="/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/mosdepth_to_download/",
    renamed_csv_path="/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv",
    output_csv_name="mosdepth_analysis_results_with_metadata.csv" # Changed output name to reflect new data
):
    """
    Performs analysis on mosdepth summary files, linking contigs to r_MAG_ids,
    counting mean > 1 occurrences, and including sample metadata and Genus.

    Args:
        summary_files_dir (str): Directory containing mosdepth *.summary.txt files.
        renamed_csv_path (str): Path to the RENAMED_CONTIG_TAXONOMY CSV file.
        output_csv_name (str): Name for the output CSV file with analysis results.
    """

    print(f"Starting analysis...")
    print(f"Summary files directory: {summary_files_dir}")
    print(f"Renamed CSV path: {renamed_csv_path}")

    # --- 1. Load RENAMED_CONTIG_TAXONOMY.csv ---
    try:
        renamed_df = pd.read_csv(renamed_csv_path)
        
        # Create mappings for quick lookup
        contig_to_mag_map = renamed_df.set_index('updated_contig_id')['r_MAG_id'].to_dict()
        
        # This will store r_MAG_id -> set of all its associated updated_contig_ids
        mag_to_all_contigs_in_sheet = {}
        for index, row in renamed_df.iterrows():
            mag_id = row['r_MAG_id']
            contig_id = row['updated_contig_id']
            if mag_id not in mag_to_all_contigs_in_sheet:
                mag_to_all_contigs_in_sheet[mag_id] = set()
            mag_to_all_contigs_in_sheet[mag_id].add(contig_id)
        
        # Create a mapping for r_MAG_id to Genus
        # Using drop_duplicates to handle cases where an r_MAG_id might appear
        # multiple times but ensuring one genus per MAG_id is picked (the first encountered)
        mag_to_genus_map = renamed_df.drop_duplicates(subset=['r_MAG_id']).set_index('r_MAG_id')['Genus'].to_dict()

        print(f"Loaded RENAMED_CONTIG_TAXONOMY.csv with {len(renamed_df)} entries.")
        print(f"Loaded {len(mag_to_all_contigs_in_sheet)} unique r_MAG_ids.")
        print(f"Loaded {len(mag_to_genus_map)} unique r_MAG_ids with Genus information.")

    except FileNotFoundError:
        print(f"Error: RENAMED_CONTIG_TAXONOMY.csv not found at {renamed_csv_path}")
        return
    except KeyError as e:
        print(f"Error: Missing expected column in RENAMED_CONTIG_TAXONOMY.csv: {e}")
        print("Please ensure 'updated_contig_id', 'r_MAG_id', and 'Genus' columns exist.")
        return
    except Exception as e:
        print(f"An unexpected error occurred loading RENAMED CSV: {e}")
        return

    # --- 2. Prepare results list ---
    analysis_results = []
    output_headers = [
        'summary_file_name',
        'State',
        'Month',
        'Caste',
        'Colony',
        'Rep',
        'r_MAG_id',
        'Genus', # Added Genus
        'contigs_observed_in_sheet',
        'contigs_observed_in_summary_txt',
        'contigs_with_mean_gt_1',
        'contigs_present_ratio'
    ]

    # --- 3. Iterate through mosdepth summary files ---
    summary_files = [f for f in os.listdir(summary_files_dir) if f.endswith('.summary.txt')]
    if not summary_files:
        print(f"No *.summary.txt files found in {summary_files_dir}. Exiting.")
        return

    print(f"\nProcessing {len(summary_files)} summary files...")

    for summary_file_name in summary_files:
        current_summary_path = os.path.join(summary_files_dir, summary_file_name)
        
        # --- Extract metadata from filename ---
        state, month, caste, colony, rep = None, None, None, None, None
        try:
            # Remove .mosdepth.summary.txt extension
            base_name = summary_file_name.replace(".mosdepth.summary.txt", "")
            parts = base_name.split('_')
            if len(parts) >= 5:
                state = parts[0]
                month = parts[1]
                caste = parts[2]
                colony = int(parts[3]) # Convert to int
                rep = int(parts[4]) # Convert to int
            else:
                print(f"    Warning: Filename '{summary_file_name}' does not match expected format (State_Month_Caste_Colony_Rep_...). Metadata will be N/A.")
        except ValueError:
            print(f"    Warning: Could not parse Colony or Rep as integer from '{summary_file_name}'. Metadata will be N/A.")
        except Exception as e:
            print(f"    Error extracting metadata from filename '{summary_file_name}': {e}. Metadata will be N/A.")

        # Store data relevant to this summary file's analysis
        # r_MAG_id -> {'contigs_in_summary': set(), 'contigs_mean_gt_1': set()}
        mags_data_for_current_file = {} 

        print(f"  Processing: {summary_file_name}")

        try:
            with open(current_summary_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader) # Skip header
                
                # Find column indices
                try:
                    chrom_idx = header.index('chrom')
                    mean_idx = header.index('mean')
                except ValueError as e:
                    print(f"    Skipping {summary_file_name}: Missing expected column in summary file header: {e}. (Expected 'chrom' and 'mean')")
                    continue

                for row in reader:
                    if len(row) > max(chrom_idx, mean_idx): # Basic check to prevent index errors
                        chrom = row[chrom_idx].strip()
                        try:
                            mean_val = float(row[mean_idx])
                        except ValueError:
                            # print(f"      Warning: Could not parse 'mean' value '{row[mean_idx]}' for chrom '{chrom}' in {summary_file_name}. Skipping row.")
                            continue # Skip row if mean is not a valid number

                        # Check if this chrom exists in our RENAMED mapping
                        if chrom in contig_to_mag_map:
                            r_mag_id = contig_to_mag_map[chrom]

                            if r_mag_id not in mags_data_for_current_file:
                                mags_data_for_current_file[r_mag_id] = {
                                    'contigs_in_summary': set(),
                                    'contigs_mean_gt_1': set()
                                }
                            
                            mags_data_for_current_file[r_mag_id]['contigs_in_summary'].add(chrom) # Track all contigs from this MAG found in this summary file
                            if mean_val > 1:
                                mags_data_for_current_file[r_mag_id]['contigs_mean_gt_1'].add(chrom)
                        else:
                            # Per instruction: "mark it but otherwise ignore it"
                            # print(f"      Warning: Contig '{chrom}' from {summary_file_name} not found in RENAMED_CONTIG_TAXONOMY.csv. Ignoring.")
                            pass # Do nothing, just move on
        except FileNotFoundError:
            print(f"    Error: Summary file not found at {current_summary_path}. Skipping.")
            continue
        except Exception as e:
            print(f"    An unexpected error occurred processing {summary_file_name}: {e}. Skipping.")
            continue

        # --- 4. Calculate metrics for each r_MAG_id found in this summary file ---
        for r_mag_id, data in mags_data_for_current_file.items():
            # Get total contigs for this MAG from the RENAMED sheet
            contigs_in_sheet = len(mag_to_all_contigs_in_sheet.get(r_mag_id, set())) # Use .get() for safety
            
            # Get Genus for this r_MAG_id
            genus = mag_to_genus_map.get(r_mag_id, "N/A") # Default to N/A if Genus not found

            contigs_in_summary_count = len(data['contigs_in_summary'])
            contigs_mean_gt_1_count = len(data['contigs_mean_gt_1'])

            contigs_present_ratio = 0.0
            if contigs_in_sheet > 0: # Avoid division by zero
                contigs_present_ratio = contigs_mean_gt_1_count / contigs_in_sheet
            
            analysis_results.append({
                'summary_file_name': summary_file_name,
                'State': state,
                'Month': month,
                'Caste': caste,
                'Colony': colony,
                'Rep': rep,
                'r_MAG_id': r_mag_id,
                'Genus': genus, # Added Genus
                'contigs_observed_in_sheet': contigs_in_sheet,
                'contigs_observed_in_summary_txt': contigs_in_summary_count,
                'contigs_with_mean_gt_1': contigs_mean_gt_1_count,
                'contigs_present_ratio': contigs_present_ratio
            })

    # --- 5. Write results to CSV ---
    if analysis_results:
        output_df = pd.DataFrame(analysis_results, columns=output_headers)
        output_path = os.path.join(os.path.dirname(summary_files_dir), output_csv_name) # Puts output one level up from mosdepth_to_download
        output_df.to_csv(output_path, index=False)
        print(f"\nAnalysis complete. Results saved to: {output_path}")
    else:
        print("\nNo analysis results to save. This might mean no relevant data was found or processed.")

if __name__ == "__main__":
    analyze_mosdepth_summaries()

    # If you want to specify different paths, uncomment and modify these lines:
    # analyze_mosdepth_summaries(
    #     summary_files_dir="/path/to/your/mosdepth_summaries/",
    #     renamed_csv_path="/path/to/your/RENAMED_CONTIG_TAXONOMY.csv",
    #     output_csv_name="my_custom_results.csv"
    # )