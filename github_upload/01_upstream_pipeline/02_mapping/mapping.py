import os
import subprocess
from pathlib import Path
import logging
from concurrent.futures import ProcessPoolExecutor

#define some paths

TRIMMED_DIR = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/01_TrimmedReads"
ASSEMBLIES_DIR = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/02_Assemblies"
MAPPING_DIR = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/03_Mapping"

#bowtie2 indexing of coassemblies
def create_index(coassembly_path, index_name):
    index_dir = coassembly_path.parent
    subprocess.run([
        "bowtie2-build",
        str(coassembly_path),
        str(index_dir / index_name)], check=True
    )
    logging.info(f"Index created for: {index_name}")
    
#actual function to map reads, process BAM files

def map_reads(hierarchy, subdir):
    #the paths and filenames
    assembly_path = Path(ASSEMBLIES_DIR) / hierarchy / "final.contigs.fa"
    index_name = f"{assembly_path.parent.name}_coassembly"
    trimmed_path = Path(TRIMMED_DIR) / hierarchy / subdir
    mapping_path = Path(MAPPING_DIR) / hierarchy / subdir
    mapping_path.mkdir(parents=True, exist_ok=True)
    
    r1_file = trimmed_path / f"{hierarchy}_{subdir}_R1.trimmed.fastq.gz"
    r2_file = trimmed_path / f"{hierarchy}_{subdir}_R2.trimmed.fastq.gz"
    sam_file = mapping_path / f"{hierarchy}_{subdir}.sam"
    bam_file = mapping_path / f"{hierarchy}_{subdir}.bam"
    sorted_bam_file = mapping_path / f"{hierarchy}_{subdir}.sorted.bam"
    
    try:
        #map reads
        logging.info(f"Mapping reads for: {hierarchy}/{subdir}")
        
        subprocess.run(
            ["bowtie2", "--very-sensitive", "-x", str(assembly_path.parent / index_name),
             "-1", str(r1_file), "-2", str(r2_file), "-S", str(sam_file)],
            check=True,
        )
        
        #sam to bam
        subprocess.run(["samtools", "view", "-bS", str(sam_file), "-o", str(bam_file)], check=True)
        
        #sort bam
        subprocess.run(["samtools", "sort", "-o", str(sorted_bam_file), str(bam_file)], check=True)
        
        #index bam
        subprocess.run(["samtools", "index", str(sorted_bam_file)], check=True)
        
        #cleanup temp files
        
        os.remove(sam_file)
        os.remove(bam_file)
        logging.info(f"Mapping completed for: {hierarchy}/{subdir}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during mapping for {hierarchy}/{subdir}: {e}")
        

def map_reads_safe(hierarchy, subdir_name):
    """
    Safely wraps the map_reads function to catch and log exceptions.
    """
    try:
        map_reads(hierarchy, subdir_name)
    except Exception as e:
        logging.error(f"Error processing {hierarchy}/{subdir_name}: {e}")

def process_directory(hierarchy):
    """
    Processes a directory by creating a Bowtie2 index and mapping reads for all subdirectories in parallel.
    """
    coassembly_path = Path(ASSEMBLIES_DIR) / hierarchy / "final.contigs.fa"
    # Use only the last part of the hierarchy as the index name
    index_name = f"{coassembly_path.parent.name}_coassembly"
    subdirs = sorted((Path(TRIMMED_DIR) / hierarchy).iterdir())

    # Create Bowtie2 index for the coassembly
    try:
        create_index(coassembly_path, index_name)
    except Exception as e:
        logging.error(f"Error creating index for {coassembly_path}: {e}")
        return

    # Process all subdirectories in parallel
    logging.info(f"Processing subdirectories for hierarchy {hierarchy}")
    try:
        with ProcessPoolExecutor(max_workers=24) as executor:
            futures = [
                executor.submit(map_reads_safe, hierarchy, subdir.name)
                for subdir in subdirs
            ]
            for future in futures:
                future.result()  # Wait for all tasks to complete and catch exceptions
    except Exception as e:
        logging.error(f"Unexpected error during parallel processing: {e}")

    logging.info(f"Finished processing {hierarchy}")
        
#main function
if __name__ == "__main__":
    hierarchies = [
        str(assembly_path.relative_to(ASSEMBLIES_DIR).parent)
        for assembly_path in Path(ASSEMBLIES_DIR).rglob("final.contigs.fa")
    ]
    for hierarchy in hierarchies:
        logging.info(f"Processing hierarchy: {hierarchy}")
        process_directory(hierarchy)
