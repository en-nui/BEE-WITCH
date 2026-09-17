import os
import subprocess

def find_and_process_bam_files(root_dir):
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(".bam"):
                bam_path = os.path.join(dirpath, filename)
                output_file = os.path.join(dirpath, f"{filename.replace('.bam', '')}_coverage.txt")
            
                if os.path.exists(output_file):
                    print(f"Skipping {bam_path}, output file already exists.")
                    continue


                command = f"bedtools genomecov -dz -ibam {bam_path} > {output_file}"
                print(f"Running: {command} on {bam_path}")
                subprocess.run(command, shell=True, check=True)

if __name__ == "__main__":
    root_directory = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/16_rMAG_mapping"
    find_and_process_bam_files(root_directory)
