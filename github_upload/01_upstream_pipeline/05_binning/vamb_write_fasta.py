import os
import subprocess
import logging

#set up log file
logging.basicConfig(
    filename="process_bins.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


def find_clusters_and_process(base_dir, assemblies_dir, create_fasta_script, size):
    #os.walk through base_dir, locate the file 'vae_clusters_unsplit.tsv', process with create_fasta.py in software/vamb/src/
    
    #args:
        #base_dir (str): base directory to look for vae_clusters_unsplit.tsv
        #assembleis_dir (str): base dir to locate 'final.contigs.fa'
        #create_fasta_script (str): path to find 'create_fasta.py'
        #size (int): size argument to pass to 'create_fasta.py' (determines size of bin to keep)

    for root, dirs, files in os.walk(base_dir):
        if "vae_clusters_unsplit.tsv" in files:
            clusterspath = os.path.join(root, "vae_clusters_unsplit.tsv")
            relative_path = os.path.relpath(root, base_dir)
            fastapath = os.path.join(assemblies_dir, relative_path, "final.contigs.fa")
            outpath = os.path.join(root, "cluster_fasta")
 
            logging.info(f"Found clusterspath: {clusterspath}")
            logging.info(f"Found fastapath: {fastapath}")
            logging.info(f"Found outpath: {outpath}")

            #if os.path.exists(outpath):
            #    logging.info(f"Outpath already exists. Skipping: {outpath}")
            #    continue

            try:
                os.makedirs(outpath, exist_ok=True)
                logging.info(f"Created outpath directory: {outpath}")
            except Exception as e:
                logging.error(f"Failed to create outpath directory {outpath}. Error: {e}")
                continue

        #run create_fasta.py
            try:
                command = [
                    "python", create_fasta_script,
                    fastapath, clusterspath,
                    str(size), outpath
                ]

                subprocess.run(command, check=True)
                logging.info(f"Successfully processed: {clusterspath}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running create_fasta.py for {clusterspath}: {e}")
            except Exception as e:
                logging.error(f"Unexpected error for {clusterspath}: {e}")



if __name__ == "__main__":
    base_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/07b_mMAG_binning_vamb"
    assemblies_dir = "/N/project/NewtonLab/robinch/projects/BEE-WITCH/02_Assemblies"
    create_fasta_script = "/N/project/NewtonLab/robinch/software/vamb/src/create_fasta.py"
    size = 250000

    find_clusters_and_process(base_dir, assemblies_dir, create_fasta_script, size)
        

