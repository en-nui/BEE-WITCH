import os
import subprocess
import glob

root_directory = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/18_anvio_outputs"

def run_variability_profile(contig_db, profile_db, splits_file, output_name, extra_flags=None):
    """Runs the anvi-gen-variability-profile command."""
    command = [
        "anvi-gen-variability-profile",
        "-p", profile_db,
        "-c", contig_db,
        "--splits-of-interest", splits_file,
        "--include-split-names",
        "--include-contig-names",
        "-o", output_name,
    ]
    if extra_flags:
        command.extend(extra_flags)

    print(f"Executing: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
        print(f"Successfully generated {output_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        return False
    except FileNotFoundError:
        print("Error: The 'anvi-gen-variability-profile' command was not found. Make sure it's in your PATH.")
        return False

# --- First pass without --engine CDN ---
print("Starting first pass of variability profile generation...")
for colony_dir_name in os.listdir(root_directory):
    colony_path = os.path.join(root_directory, colony_dir_name)
    if os.path.isdir(colony_path):
        print(f"\nProcessing directory: {colony_path}")

        # Check for checkpoint file
        checkpoint_file = os.path.join(colony_path, "variability_done.txt")
        if os.path.exists(checkpoint_file):
            print(f"Checkpoint found for {colony_path}. Skipping first pass.")
            continue

        # Find contig database
        contig_dbs = glob.glob(os.path.join(colony_path, "*concat.db"))
        if not contig_dbs:
            print(f"Error: No *concat.db file found in {colony_path}")
            continue
        contig_db = contig_dbs[0]

        # Find profile database
        profile_subdirs = glob.glob(os.path.join(colony_path, "*_merged-profile"))
        if not profile_subdirs:
            print(f"Error: No *_merged-profile directory found in {colony_path}")
            continue
        profile_db = os.path.join(profile_subdirs[0], "PROFILE.db")
        if not os.path.exists(profile_db):
            print(f"Error: PROFILE.db not found in {profile_subdirs[0]}")
            continue

        # Find splits file
        splits_file = os.path.join(colony_path, "splits_to_import.txt")
        if not os.path.exists(splits_file):
            print(f"Error: splits_to_import.txt not found in {colony_path}")
            continue

        # Define output name
        output_name = os.path.join(colony_path, f"{colony_dir_name}_variability.txt")

        # Run the command
        if run_variability_profile(contig_db, profile_db, splits_file, output_name):
            # Create checkpoint file
            with open(checkpoint_file, "w") as f:
                f.write("Done")
            print(f"Checkpoint created: {checkpoint_file}")

print("\nFirst pass complete.")

# --- Second pass with --engine CDN ---
print("\nStarting second pass of variability profile generation with --engine CDN...")
for colony_dir_name in os.listdir(root_directory):
    colony_path = os.path.join(root_directory, colony_dir_name)
    if os.path.isdir(colony_path):
        print(f"\nProcessing directory: {colony_path}")

        # Check for first pass checkpoint file
        checkpoint_file_pass1 = os.path.join(colony_path, "variability_done.txt")
        if not os.path.exists(checkpoint_file_pass1):
            print(f"First pass not completed for {colony_path}. Skipping second pass.")
            continue

        # Check for second pass checkpoint file
        checkpoint_file_pass2 = os.path.join(colony_path, "variability_cdn_done.txt")
        if os.path.exists(checkpoint_file_pass2):
            print(f"Checkpoint found for CDN pass in {colony_path}. Skipping second pass.")
            continue

        # Find contig database (assuming it's the same as before)
        contig_dbs = glob.glob(os.path.join(colony_path, "*concat.db"))
        if not contig_dbs:
            print(f"Error: No *concat.db file found in {colony_path}")
            continue
        contig_db = contig_dbs[0]

        # Find profile database (assuming it's the same as before)
        profile_subdirs = glob.glob(os.path.join(colony_path, "*_merged-profile"))
        if not profile_subdirs:
            print(f"Error: No *_merged-profile directory found in {colony_path}")
            continue
        profile_db = os.path.join(profile_subdirs[0], "PROFILE.db")
        if not os.path.exists(profile_db):
            print(f"Error: PROFILE.db not found in {profile_subdirs[0]}")
            continue

        # Find splits file (assuming it's the same as before)
        splits_file = os.path.join(colony_path, "splits_to_import.txt")
        if not os.path.exists(splits_file):
            print(f"Error: splits_to_import.txt not found in {colony_path}")
            continue

        # Define output name for CDN run
        output_name_cdn = os.path.join(colony_path, f"{colony_dir_name}_variability_CDN.txt")

        # Run the command with --engine CDN
        if run_variability_profile(contig_db, profile_db, splits_file, output_name_cdn, extra_flags=["--engine", "CDN"]):
            # Create checkpoint file for CDN run
            with open(checkpoint_file_pass2, "w") as f:
                f.write("Done")
            print(f"Checkpoint created: {checkpoint_file_pass2}")

print("\nSecond pass with --engine CDN complete.")
