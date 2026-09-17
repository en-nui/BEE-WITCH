import pandas as pd
import os
import re

def parse_pfam2go(mapping_file):
    """
    Parses the pfam2go.txt file into a clean pandas DataFrame.
    """
    print(f"Parsing the GO mapping file: {mapping_file}")
    
    # Regex to capture the four desired components from a line
    # 1: Pfam ID (e.g., PF00001)
    # 2: Pfam Label (e.g., 7tm_1)
    # 3: GO Label (e.g., G protein-coupled receptor signaling pathway)
    # 4: GO ID (e.g., 0007186)
    line_regex = re.compile(r"Pfam:(PF\d{5})\s(.*?)\s>\sGO:(.*?)\s;\sGO:(\d+)")
    
    parsed_data = []
    try:
        with open(mapping_file, 'r') as f:
            for line in f:
                match = line_regex.search(line)
                if match:
                    pfam_id, pfam_label, go_label, go_id = match.groups()
                    parsed_data.append({
                        'pfam_id': pfam_id.strip(),
                        'pfam_label': pfam_label.strip(),
                        'GO_label': go_label.strip(),
                        'GO_id': go_id.strip()
                    })
    except FileNotFoundError:
        print(f"Error: The mapping file {mapping_file} was not found.")
        return None

    if not parsed_data:
        print("Warning: No data was successfully parsed from the mapping file.")
        return None
        
    # Since one Pfam ID can have multiple GO terms, we group by Pfam ID
    # and aggregate the GO terms into lists. This prevents row duplication during the merge.
    go_map_df = pd.DataFrame(parsed_data)
    
    # Aggregate GO terms into lists for each Pfam ID
    agg_funcs = {
        'pfam_label': 'first', # pfam_label is the same for each pfam_id
        'GO_label': lambda x: list(x),
        'GO_id': lambda x: list(x)
    }
    go_map_aggregated = go_map_df.groupby('pfam_id').agg(agg_funcs).reset_index()

    print(f"Successfully created a lookup map with {len(go_map_aggregated)} unique Pfam IDs.")
    return go_map_aggregated

def annotate_files_with_go(input_dir, output_dir, go_map):
    """
    Loops through data files, extracts Pfam IDs, and merges GO annotations.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    for filename in os.listdir(input_dir):
        if filename.endswith(".csv"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename.replace(".csv", "_annotated.csv"))
            
            print(f"Processing file: {filename}")
            
            # Load the data file
            data_df = pd.read_csv(input_path)
            
            # Extract the core Pfam ID (e.g., 'PF13936') from the 'gene_name' column
            # This regex looks for '[PF' followed by 5 digits and captures them.
            data_df['pfam_id'] = data_df['gene_name'].str.extract(r'\[(PF\d{5})')
            
            # Merge the GO annotations using the new pfam_id column
            annotated_df = pd.merge(data_df, go_map, on='pfam_id', how='left')
            
            # Save the new, enriched dataframe
            annotated_df.to_csv(output_path, index=False)
            print(f"  -> Saved annotated file to: {output_path}")

# --- HOW TO USE ---
if __name__ == '__main__':
    # 1. Define the path to your pfam2go mapping file
    pfam_map_file = "/home/robinch/projects/BEE-WITCH/pfam2go_dir/pfam2go.txt"
    
    # 2. Define the directory containing your data files to be annotated
    data_directory = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/annotated_with_unique_gene_id"
    
    # 3. Define where to save the new annotated files
    output_directory = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/final_annotated_files"
    
    # Run the workflow
    go_lookup_table = parse_pfam2go(pfam_map_file)
    if go_lookup_table is not None:
        annotate_files_with_go(data_directory, output_directory, go_lookup_table)
        print("\n--- Annotation process complete! ---")