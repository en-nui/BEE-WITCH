import pandas as pd
import os
import re

# --- Part 1: Updated Parsing Functions ---

def parse_go_slim(slim_file_path):
    """Parses a GO slim OBO file to extract IDs, names, and namespaces."""
    print(f"Parsing GO slim file: {slim_file_path}")
    slim_ids = set()
    name_map = {}
    namespace_map = {}
    current_id = None
    try:
        with open(slim_file_path, 'r') as f:
            for line in f:
                if line.startswith("id: GO:"):
                    current_id = line.strip().split(" ")[1]
                    slim_ids.add(current_id)
                elif line.startswith("name:") and current_id:
                    name_map[current_id] = line.strip().split("name: ")[1]
                elif line.startswith("namespace:") and current_id:
                    namespace_map[current_id] = line.strip().split("namespace: ")[1]
    except FileNotFoundError:
        print(f"Error: The file {slim_file_path} was not found.")
        return None, None, None
    print(f"  -> Found {len(slim_ids)} terms in the prokaryote GO slim.")
    return slim_ids, name_map, namespace_map

def parse_full_go_ontology(obo_file_path):
    """Parses go-basic.obo for an ancestry map, name map, and namespace map."""
    print(f"Parsing full GO ontology file: {obo_file_path}")
    ancestry_map = {}
    name_map = {}
    namespace_map = {}
    current_term_id = None
    try:
        with open(obo_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "[Term]":
                    current_term_id = None
                elif line.startswith("id: GO:"):
                    current_term_id = line.split(" ")[1]
                    if current_term_id not in ancestry_map:
                        ancestry_map[current_term_id] = []
                elif line.startswith("name:") and current_term_id:
                    name_map[current_term_id] = line.split("name: ")[1]
                elif line.startswith("namespace:") and current_term_id:
                    namespace_map[current_term_id] = line.split("namespace: ")[1]
                elif line.startswith("is_a: GO:") and current_term_id:
                    parent_id = line.split(" ")[1]
                    ancestry_map[current_term_id].append(parent_id)
    except FileNotFoundError:
        print(f"Error: The file {obo_file_path} was not found.")
        return None, None, None
    print(f"  -> Built maps for {len(ancestry_map)} GO terms.")
    return ancestry_map, name_map, namespace_map

# --- Part 2: Ancestry Tracing Function (unchanged) ---

def find_slim_ancestor(go_id, ancestry_map, slim_term_set, cache):
    # ... (This function is the same as the previous version)
    if go_id in cache: return cache[go_id]
    if go_id in slim_term_set: cache[go_id] = go_id; return go_id
    if go_id not in ancestry_map: cache[go_id] = None; return None
    for parent_id in ancestry_map[go_id]:
        result = find_slim_ancestor(parent_id, ancestry_map, slim_term_set, cache)
        if result: cache[go_id] = result; return result
    cache[go_id] = None; return None

# --- Part 3: Main Workflow ---

if __name__ == '__main__':
    # a. Define file paths
    full_obo_file = "/home/robinch/databases/go_basic/go-basic.obo"
    slim_obo_file = "/home/robinch/databases/go_basic/goslim_prokaryote.obo"
    annotated_files_dir = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/final_annotated_files"
    output_dir = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/slim_annotated_files"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # b. Run parsing functions
    slim_term_set, slim_name_map, slim_namespace_map = parse_go_slim(slim_obo_file)
    ancestry_map, go_name_map, go_namespace_map = parse_full_go_ontology(full_obo_file)
    
    if slim_term_set and ancestry_map:
        # c. Gather all unique GO IDs from data files
        print("\nGathering all unique GO IDs...")
        all_go_ids = set()
        for filename in os.listdir(annotated_files_dir):
            if filename.endswith("_annotated.csv"):
                filepath = os.path.join(annotated_files_dir, filename)
                df = pd.read_csv(filepath, low_memory=False)
                df_subset = df.dropna(subset=['GO_id'])
                if not df_subset.empty:
                    exploded_ids = df_subset['GO_id'].apply(eval).explode()
                    standardized_ids = "GO:" + exploded_ids.astype(str).str.zfill(7)
                    all_go_ids.update(standardized_ids.unique())
        print(f"  -> Found {len(all_go_ids)} unique GO IDs.")
        
        # d. Build the specific-to-slim map
        print("\nBuilding the specific-to-slim GO map...")
        specific_to_slim_map = {}
        memoization_cache = {}
        for go_id in all_go_ids:
            slim_ancestor = find_slim_ancestor(go_id, ancestry_map, slim_term_set, memoization_cache)
            if slim_ancestor:
                specific_to_slim_map[go_id] = slim_ancestor
        print(f"  -> Successfully mapped {len(specific_to_slim_map)} terms.")
        
        # e. Create a final, detailed mapping DataFrame for merging
        print("\nCreating detailed mapping table...")
        map_data = []
        for specific_id, slim_id in specific_to_slim_map.items():
            map_data.append({
                'GO_id_specific': specific_id,
                'GO_slim_id': slim_id,
                'GO_slim_name': slim_name_map.get(slim_id),
                'GO_slim_namespace': slim_namespace_map.get(slim_id)
            })
        final_mapping_df = pd.DataFrame(map_data)

        # f. Loop through files again to merge and save final output
        print("\nMerging slim annotations into final files...")
        for filename in os.listdir(annotated_files_dir):
            if filename.endswith("_annotated.csv"):
                filepath = os.path.join(annotated_files_dir, filename)
                df = pd.read_csv(filepath, low_memory=False)

                # Explode and standardize GO_id to prepare for merge
                df_subset = df.dropna(subset=['GO_id'])
                exploded = df_subset.assign(GO_id_specific_list=df_subset['GO_id'].apply(eval)).explode('GO_id_specific_list')
                exploded['GO_id_specific'] = "GO:" + exploded['GO_id_specific_list'].astype(str).str.zfill(7)
                
                # Merge with the slim mapping table
                merged_df = pd.merge(exploded, final_mapping_df, on='GO_id_specific', how='left')

                # Save the final, fully enriched file
                output_path = os.path.join(output_dir, filename.replace("_annotated.csv", "_slim_annotated.csv"))
                merged_df.to_csv(output_path, index=False)
                print(f"  -> Saved {output_path}")

        print("\n--- Final annotation process complete! ---")