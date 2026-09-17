import pandas as pd
import os
import glob
import logging
import time

# --- Configuration ---
# Base directories for input and output files.
BASE_ANVIO_DIR = '/home/robinch/projects/BEE-WITCH/22_anvio'
MAPPING_OUTPUT_DIR = '/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output'
OUTPUT_DIR = '/home/robinch/projects/BEE-WITCH/34_mMAG_popgen'

# --- Helper Functions ---

def setup_logger(log_file, log_to_console=True):
    """
    Configures a logger to write to a specific file and optionally to the console.
    """
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    # Clear existing handlers to prevent duplicate log messages
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        
    handlers = [logging.FileHandler(log_file)]
    if log_to_console:
        handlers.append(logging.StreamHandler())
        
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger()

def create_filtering_lookup_table_for_colony(variability_file, logger):
    """
    Creates a set of unique (sample_id, contig_name) tuples for a single colony
    that pass the filtering criteria. This is done once per colony for efficiency.
    """
    logger.info(f"Creating filtering lookup table for {os.path.basename(variability_file)}...")
    
    try:
        # Load all unique (sample_id, contig_name) pairs from the variability file.
        # This is the most memory-intensive step, but it's much smaller than the full file.
        variability_unique_df = pd.read_csv(
            variability_file, 
            sep='\t', 
            usecols=['sample_id', 'contig_name']
        ).drop_duplicates()
        
        logger.info(f"Found {len(variability_unique_df)} unique (sample_id, contig_name) pairs to filter.")
    except Exception as e:
        logger.error(f"Error reading unique pairs from {variability_file}: {e}")
        return None

    # This set will hold the final, deduplicated, and validated pairs
    valid_contig_pairs = set()

    # Iterate over the unique pairs DataFrame
    for sample_id, group in variability_unique_df.groupby('sample_id'):
        mmag_file = os.path.join(MAPPING_OUTPUT_DIR, f'{sample_id}.filtered_mmag_contigs.csv')
        
        if not os.path.exists(mmag_file):
            logger.warning(f"mMAG file not found for sample_id '{sample_id}'. Skipping filtering for this sample.")
            continue
            
        try:
            # Load the updated_contig_id list for this specific sample
            mmag_df = pd.read_csv(mmag_file, usecols=['updated_contig_id'])
            updated_contig_ids = set(mmag_df['updated_contig_id'])
            
            # Filter the current group of contigs against the mMAG list
            matching_contigs = group[group['contig_name'].isin(updated_contig_ids)]
            
            # Add the valid pairs to our final set. The set ensures no duplicates.
            for _, row in matching_contigs.iterrows():
                valid_contig_pairs.add((row['sample_id'], row['contig_name']))
                
            logger.info(f"Processed sample '{sample_id}'. Found {len(matching_contigs)} valid contigs.")
            
        except Exception as e:
            logger.error(f"Error processing mMAG file '{mmag_file}': {e}. Skipping this sample_id.")
            continue
    
    logger.info(f"Lookup table for this colony contains {len(valid_contig_pairs)} valid pairs.")
    return valid_contig_pairs

def process_colony(colony_path, logger):
    """
    Main processing logic for a single colony. It reads the large variability
    file line by line and filters based on the pre-built lookup table for this colony.
    """
    colony_id = os.path.basename(colony_path).split('_')[1]
    logger.info(f"--- Starting processing for Colony {colony_id} ---")

    # 1. Define file paths - this is now dynamic
    variability_file = os.path.join(colony_path, f'{colony_id}_merged_variability.txt')
    colony_output_dir = os.path.join(OUTPUT_DIR, f'colony_{colony_id}')
    output_path = os.path.join(colony_output_dir, f'filtered_variability.txt')
    
    if not os.path.exists(variability_file):
        logger.error(f"Variability file not found for Colony {colony_id}. Skipping.")
        return
        
    os.makedirs(colony_output_dir, exist_ok=True)
    
    # Create the filtering lookup table specifically for THIS colony
    valid_contig_pairs = create_filtering_lookup_table_for_colony(variability_file, logger)
    if not valid_contig_pairs:
        logger.error("Failed to create filtering lookup table for this colony. Exiting.")
        return

    # 2. Iterate through the variability file and write filtered data to output
    logger.info("Iterating through the full variability file and writing filtered data...")
    
    total_rows = 0
    filtered_rows = 0
    
    try:
        with open(variability_file, 'r') as infile, open(output_path, 'w') as outfile:
            # Read header and write it to the output file
            header = infile.readline().strip().split('\t')
            outfile.write('\t'.join(header) + '\n')
            
            # Process remaining lines
            for line in infile:
                total_rows += 1
                parts = line.strip().split('\t')
                
                # Use column names to find indices
                try:
                    sample_id_index = header.index('sample_id')
                    contig_name_index = header.index('contig_name')
                except ValueError:
                    logger.error("Required columns 'sample_id' or 'contig_name' not found in header. Skipping file.")
                    break
                
                # Check if we have enough parts and if the lookup is valid
                if len(parts) > max(sample_id_index, contig_name_index):
                    sample_id = parts[sample_id_index]
                    contig_name = parts[contig_name_index]
                    
                    # Fast lookup against our pre-built set
                    if (sample_id, contig_name) in valid_contig_pairs:
                        outfile.write(line)
                        filtered_rows += 1
                
                if total_rows % 1000000 == 0:
                    logger.info(f"Processed {total_rows} rows. Found {filtered_rows} matching rows.")
                    
    except Exception as e:
        logger.error(f"Error while streaming variability file: {e}. Output may be incomplete.")
        return
        
    logger.info(f"Finished processing. Total rows read: {total_rows}. Total filtered rows: {filtered_rows}.")
    logger.info(f"Successfully saved filtered variability data to {output_path}")
    logger.info(f"--- Finished processing for Colony {colony_id} ---\n")


def main():
    """
    Main function to orchestrate the entire workflow.
    """
    start_time = time.time()
    
    log_file_path = os.path.join(OUTPUT_DIR, 'mmag_filtering.log')
    main_logger = setup_logger(log_file_path)
    main_logger.info("Starting the mMAG filtering workflow...")
    
    colony_dirs = glob.glob(os.path.join(BASE_ANVIO_DIR, 'Colony_*'))
    
    if not colony_dirs:
        main_logger.error(f"No 'Colony_*' directories found in {BASE_ANVIO_DIR}.")
        return

    main_logger.info(f"Found {len(colony_dirs)} colony directories to process.")

    for colony_path in sorted(colony_dirs):
        process_colony(colony_path, main_logger)
        
    end_time = time.time()
    main_logger.info(f"Workflow finished. Total time: {end_time - start_time:.2f} seconds.")


if __name__ == '__main__':
    main()

