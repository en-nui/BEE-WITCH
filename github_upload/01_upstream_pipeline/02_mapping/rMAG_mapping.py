import os
import glob
import subprocess
import logging
import multiprocessing

# Configure logging
log_file = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/16_rMAG_mapping/mapping_log.txt"
logging.basicConfig(filename=log_file, level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(cmd):
    try:
        process = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"Command executed successfully: {cmd}")
        return process.stdout.decode('utf-8')
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        logging.error(f"Return code: {e.returncode}")
        logging.error(f"Stdout: {e.stdout.decode('utf-8')}")
        logging.error(f"Stderr: {e.stderr.decode('utf-8')}")
        return None
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        return None

def find_read_dirs(colony_id):
    matching_dirs = glob.glob(f"{reads_base_dir}/*/*/*/{colony_id}/*", recursive=True)
    return matching_dirs

def index_fna(fna_file, bowtie2_threads): # Pass bowtie2_threads as argument
    index_name = os.path.splitext(os.path.basename(fna_file))[0]
    cmd = f"bowtie2-build {fna_file} {concat_fna_dir}/{index_name}" # Use bowtie2_threads
    run_command(cmd)
    logging.info(f"Indexed: {fna_file} -> {index_name}")

def map_reads(fna_file, bowtie2_threads, samtools_threads): # Pass thread counts as arguments
    index_name = os.path.splitext(os.path.basename(fna_file))[0]
    colony_id = index_name.replace("_concat", "")
    read_dirs = find_read_dirs(colony_id)

    for read_dir in read_dirs:
        output_dir = read_dir.replace(reads_base_dir, output_base_dir)
        os.makedirs(output_dir, exist_ok=True)

        r1_files = glob.glob(f"{read_dir}/*_R1.trimmed.fastq.gz")
        r2_files = glob.glob(f"{read_dir}/*_R2.trimmed.fastq.gz")

        if not r1_files or not r2_files:
            continue

        paired_reads = []
        for r1 in r1_files:
            base_name = r1.replace('R1.trimmed.fastq.gz', '')
            r2 = f"{base_name}R2.trimmed.fastq.gz"
            if r2 in r2_files:
                paired_reads.append((r1, r2))

        if not paired_reads:
            continue

        r1, r2 = paired_reads[0]
        sam_base_name = os.path.basename(r1).replace('_R1.trimmed.fastq.gz', '')
        sam_file = os.path.join(output_dir, f"{sam_base_name}.sam")
        bam_file = os.path.join(output_dir, f"{sam_base_name}.bam")

        cmd_map = f"bowtie2 -p {bowtie2_threads} -x {concat_fna_dir}/{index_name} --very-sensitive -1 {r1} -2 {r2} -S {sam_file}" # Use bowtie2_threads
        run_command(cmd_map)

        cmd_bam = f"samtools view -bS {sam_file} > {bam_file}"
        bam_result = run_command(cmd_bam)

        if bam_result is None:
            logging.error(f"Failed to create BAM: {bam_file}")
            continue

        cmd_sort = f"samtools sort -@ {samtools_threads} {bam_file} -o {bam_file}.sorted" # Use samtools_threads
        sort_result = run_command(cmd_sort)

        if sort_result is None:
            logging.error(f"Failed to sort BAM: {bam_file}")
            continue

        os.rename(f"{bam_file}.sorted", bam_file)

        cmd_index = f"samtools index {bam_file}"
        run_command(cmd_index)

        os.remove(sam_file)
        logging.info(f"Mapped: {r1}, {r2} -> {bam_file}")


def process_colony(fna_file, bowtie2_threads, samtools_threads): # Modified process_colony to accept thread counts
    index_fna(fna_file, bowtie2_threads) # Pass bowtie2_threads
    map_reads(fna_file, bowtie2_threads, samtools_threads) # Pass thread counts

if __name__ == "__main__":
    concat_fna_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/15_MAG_dRep/colony_concatenations"
    reads_base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/01_TrimmedReads"
    output_base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/16_rMAG_mapping"
    os.makedirs(output_base_dir, exist_ok=True)

    fna_files = glob.glob(os.path.join(concat_fna_dir, "*_concat.fna"))

    num_processes = 4 # Match --ntasks-per-node
    bowtie2_threads_per_process = 2 # Match --cpus-per-task
    samtools_threads_per_process = 2 # Match --cpus-per-task

    logging.info(f"Using {num_processes} processes for parallel mapping.")
    logging.info(f"bowtie2 threads per process: {bowtie2_threads_per_process}")
    logging.info(f"samtools sort threads per process: {samtools_threads_per_process}")


    with multiprocessing.Pool(processes=num_processes) as pool:
        # Pass thread counts to process_colony
        pool.starmap(process_colony, [(fna_file, bowtie2_threads_per_process, samtools_threads_per_process) for fna_file in fna_files])

    logging.info("Parallel mapping completed for all colonies.")
    print("Mapping process initiated for all colonies in parallel. Check the log file for progress.")
