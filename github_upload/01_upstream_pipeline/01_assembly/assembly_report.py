#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

# Base directories
BASE_TRIMMED_READS = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/01_TrimmedReads"
BASE_ASSEMBLIES = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/03_Assemblies"
OUTPUT_DIR = "./report_outputs"  # Directory to store the report files

def check_assembly_status(trimmed_reads_base, assemblies_base, output_dir):
    # Output files
    status_report_file = os.path.join(output_dir, "assembly_status_report.txt")
    failed_assemblies_file = os.path.join(output_dir, "failed_assemblies.txt")
    not_started_assemblies_file = os.path.join(output_dir, "not_started_assemblies.txt")
    completed_assemblies_file = os.path.join(output_dir, "completed_assemblies.txt")

    # Initialize outputs
    status_report = []
    failed_assemblies = []
    not_started_assemblies = []
    completed_assemblies = []

    # Traverse the trimmed reads directory at the coassembly level
    for dirpath, dirnames, filenames in os.walk(trimmed_reads_base):
        # Check if temp_R1.fq.gz and temp_R2.fq.gz are present in this directory
        temp_files = {"temp_R1.fq.gz", "temp_R2.fq.gz"}
        found_temp_files = temp_files.issubset(set(filenames))

        # Derive the corresponding assembly directory
        relative_path = os.path.relpath(dirpath, trimmed_reads_base)
        assembly_dir = os.path.join(assemblies_base, relative_path)
        depth = len(relative_path.split(os.sep))
    
        if depth != 4:
            continue

        if found_temp_files:
            # Check if the assembly directory exists
            if os.path.exists(assembly_dir):
                final_contigs_file = os.path.join(assembly_dir, "final.contigs.fa")
                intermediate_contigs_dir = os.path.join(assembly_dir, "intermediate_contigs")

                if os.path.isfile(final_contigs_file):
                    # Assembly is complete
                    status_report.append(f"{relative_path}: Completed")
                    completed_assemblies.append(relative_path)
                elif os.path.isdir(intermediate_contigs_dir) and os.listdir(intermediate_contigs_dir):
                    # Intermediate contigs exist, but final contigs are missing
                    status_report.append(f"{relative_path}: In Progress")
                else:
                    # Intermediate contigs directory is empty
                    status_report.append(f"{relative_path}: Failed")
                    failed_assemblies.append(relative_path)
            else:
                # Assembly has not started
                status_report.append(f"{relative_path}: Not Yet Started")
                not_started_assemblies.append(relative_path)
        else:
            # No temp files found, assembly has not started
            status_report.append(f"{relative_path}: Not Yet Started")
            not_started_assemblies.append(relative_path)

    # Write the results to output files
    with open(status_report_file, "w") as f:
        f.write("\n".join(status_report))

    with open(failed_assemblies_file, "w") as f:
        f.write("\n".join(failed_assemblies))

    with open(not_started_assemblies_file, "w") as f:
        f.write("\n".join(not_started_assemblies))
    
    with open(completed_assemblies_file, "w") as f:
        f.write("\n".join(completed_assemblies))

    print(f"Status report saved to {status_report_file}")
    print(f"Failed assemblies saved to {failed_assemblies_file}")
    print(f"Completed assemblies saved to {completed_assemblies_file}")
    print(f"Not started assemblies saved to {not_started_assemblies_file}")


if __name__ == "__main__":
    trimmed_reads_base = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/01_TrimmedReads"
    assemblies_base = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/03_Assemblies"
    output_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/misc/assembly_progress_reports"

    os.makedirs(output_dir, exist_ok=True)

    check_assembly_status(trimmed_reads_base, assemblies_base, output_dir)

