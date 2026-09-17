import os
import subprocess
import logging
from concurrent.futures import ProcessPoolExecutor


# Directories
raw_reads_dir = '/N/project/NewtonLab/robinch/projects/BEE-WITCH/00_RawReads'
trimmed_reads_dir = '01_TrimmedReads'

# Set up logging
logging.basicConfig(level=logging.INFO)

def trim_and_quality_score(paired_files):
    """
    Trim adapters and quality-score filter a FastQ file using BBDuk.
    Paired files (R1 and R2) are processed together.
    """
    r1_file, r2_file = paired_files  # Unpack the paired files tuple

    try:
        # Create output directory for trimmed files, preserving the directory structure
        relative_path = os.path.relpath(os.path.dirname(r1_file), raw_reads_dir)
        output_dir = os.path.join(trimmed_reads_dir, relative_path)
        os.makedirs(output_dir, exist_ok=True)

        # Define the output filenames for trimmed R1 and R2 files
        output_r1 = os.path.join(output_dir, os.path.basename(r1_file).replace('.fastq.gz', '.trimmed.fastq.gz'))
        output_r2 = os.path.join(output_dir, os.path.basename(r2_file).replace('.fastq.gz', '.trimmed.fastq.gz'))

        logging.info(f"Processing paired files: {r1_file} and {r2_file}")

        # Run bbduk.sh for paired-end trimming and quality filtering
        subprocess.run(
            ['bbduk.sh', f'in1={r1_file}', f'in2={r2_file}', f'out1={output_r1}', f'out2={output_r2}',  
             'ref=/N/soft/rhel7/bbtools/38.72/resources/adapters.fa', 'ktrim=r', 'k=23', 'mink=11', 
             'hdist=1', 'minlen=50', 'tpe', 'tbo', 'qtrim=r', 'trimq=10', 'ftm=5', 'ftl=10'],
            check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )
        logging.info(f"Trimmed and quality scored: {r1_file} -> {output_r1}, {r2_file} -> {output_r2} ")
    except Exception as e:
        logging.error(f"Error trimming files {r1_file} and {r2_file}: {e}")


def process_fastq_files():
    """
    Trim and quality-score all FastQ files in parallel.
    """
    paired_files = []  # List to hold tuples of (R1_file, R2_file)
    
    for root, _, files in os.walk(raw_reads_dir):
        r1_files = [file for file in files if file.endswith('_R1.fastq.gz')]
        for r1_file in r1_files:
            r2_file = r1_file.replace('_R1.fastq.gz', '_R2.fastq.gz')
            if r2_file in files:
                paired_files.append((os.path.join(root, r1_file), os.path.join(root, r2_file)))
    
    with ProcessPoolExecutor(max_workers=24) as executor:
        executor.map(trim_and_quality_score, paired_files)

    logging.info("All FastQ files processed successfully.")
if __name__ == '__main__':
    process_fastq_files()

