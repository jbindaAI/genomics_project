from concurrent.futures import ProcessPoolExecutor
import subprocess
import argparse
from pathlib import Path
from tqdm import tqdm
import os
import re


def add_paralog_identifiers(msa_file_path: str, new_dir:str) -> str:
    """
    Parses the MSA file to add a number to each duplicated genome name.
    """
    # Create a dictionary to keep track of duplicate names
    name_counts = {}
    output_lines = []
    
    with open(msa_file_path, 'r') as msa_file:
        for line in msa_file:
            if line.startswith(">"):  # Fasta header lines begin with '>'
                genome_name = line.strip().lstrip(">")  # Extract genome name
                if genome_name in name_counts:
                    name_counts[genome_name] += 1
                    new_name = f"{genome_name}_{name_counts[genome_name]}"
                else:
                    name_counts[genome_name] = 0
                    new_name = genome_name
                # Replace the genome name in the header line
                output_lines.append(f">{new_name}\n")
            else:
                # Preserve the sequence lines as they are
                output_lines.append(line)
    
    # Save the modified MSA to a new file
    new_msa_file_name = msa_file_path.split("/")[-1].replace(".fasta", "_modified.msa")
    new_msa_file_path = os.path.join(new_dir, new_msa_file_name)
    with open(new_msa_file_path, 'w') as new_msa_file:
        new_msa_file.writelines(output_lines)
    return new_msa_file_path.split("/")[-1]


def run_tree_computation(msa_file: str, msa_path: str, output_dir: str, cpu_cores: int, bootstrap: int):
    exact_filepath = os.path.join(msa_path, msa_file)
    output_prefix = os.path.join(output_dir, f"{msa_file.split('-')[0]}")
    
    try:
        # Prepare the IQ-TREE command
        command = [
            "iqtree",
            "-s", exact_filepath,
            "-T", str(cpu_cores),
            "-pre", output_prefix,
            "-m", "WAG+I",
            "-quiet"
        ]

        if bootstrap > 0:
            command.extend(["-b", str(bootstrap)])
        
        # Run the IQ-TREE command
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Return success
        return (msa_file, "Success")
    except subprocess.CalledProcessError as e:
        return (msa_file, f"Failed (Code {e.returncode})")
    except FileNotFoundError as e:
        return (msa_file, "File Not Found")
    except Exception as e:
        return (msa_file, f"Unexpected Error: {e}")


def make_trees_paralogs(msa_path: str, output_dir: str, cpu_cores: int, bootstrap: int, num_processes: int = 4):
    os.makedirs(output_dir, exist_ok=True)
    mod_msa_path = "allignment/msa_results/paralogs_modified/bacteria/"
    os.makedirs(mod_msa_path, exist_ok=True) 

    msa_files = os.listdir(msa_path)
    msa_files = [add_paralog_identifiers(os.path.join(msa_path, msa_file), mod_msa_path) for msa_file in msa_files]
    
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        results = list(tqdm(
            executor.map(
                run_tree_computation,
                msa_files,
                [mod_msa_path] * len(msa_files),
                [output_dir] * len(msa_files),
                [cpu_cores // num_processes] * len(msa_files),
                [bootstrap] * len(msa_files)
            ),
            total=len(msa_files),
            desc="Computing trees",
            unit="file"
        ))


def make_trees_ortologs(msa_path: str, output_dir: str, cpu_cores: int, bootstrap: int, num_processes: int = 4):
    os.makedirs(output_dir, exist_ok=True)
    
    msa_files = os.listdir(msa_path)
    
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        results = list(tqdm(
            executor.map(
                run_tree_computation,
                msa_files,
                [msa_path] * len(msa_files),
                [output_dir] * len(msa_files),
                [cpu_cores // num_processes] * len(msa_files),
                [bootstrap] * len(msa_files)
            ),
            total=len(msa_files),
            desc="Computing trees",
            unit="file"
        ))


def compute_bootstrap_support(tree: str) -> float:
    """
    Computes the average bootstrap support for a given tree.
    """
    bootstrap_values = re.findall(r'\)(\d+):', tree)
    bootstrap_values = [float(val) for val in bootstrap_values]

    # Compute the average bootstrap support
    if bootstrap_values:
        average_support = sum(bootstrap_values) / len(bootstrap_values)
    else:
        average_support = 0.0
    
    return average_support


def merge_results(trees_path: str, output_file: str = "all_trees.txt", eliminate_trees:bool=False, support_threshold:float=70) -> None:
    """
    Merges all .treefile files in the given directory into one file.
    """
    trees_dir = Path(trees_path)
    treefiles = sorted(trees_dir.glob("*.treefile"))

    output_path = trees_dir / output_file
    try:
        with output_path.open("w") as f_out:
            for treefile in treefiles:
                with treefile.open("r") as f_in:
                    first_line = f_in.readline().strip() # tree
                    if first_line:
                        if eliminate_trees:
                            support = compute_bootstrap_support(first_line)
                            if support < support_threshold:
                                continue
                        f_out.write(first_line + "\n")
        print(f"Merged {len(treefiles)} .treefile(s) into {output_file}")
    except Exception as e:
        print(f"Error while merging files: {e}")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for performing tree computations.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    parser.add_argument("--cpu_cores", required=True, type=int, help="Number of CPU cores to use during computations.")
    parser.add_argument("--num_processes", required=True, type=int, help="Number of separate processes to run. CPU cores are uniformly distributed on processes.")
    parser.add_argument("--bootstrap", required=True, type=int, help="Number of bootstrap replicates.")
    parser.add_argument("--support_threshold", required=True, type=float, help="Threshold for mean bootstrap support in Tree.")
    args = parser.parse_args()

    CPU_CORES = args.cpu_cores
    NUM_PROCESSES = args.num_processes
    BOOTSTRAP = args.bootstrap
    SUPPORT_THRESHOLD = args.support_threshold
    BASENAME = args.basename

    ORTOLOGS_MSA_PATH = os.path.join("allignment/msa_results/ortologs", BASENAME)
    PARALOGS_MSA_RESULTS = os.path.join("allignment/msa_results/paralogs", BASENAME)

    ORTOLOGS_OUTPUT_DIR = os.path.join("trees/tree_results/", BASENAME, "ortologs")
    PARALOGS_OUTPUT_DIR = os.path.join("trees/tree_results/", BASENAME, "paralogs")

    ORTOLOGS_BOOTSTRAP_OUTPUT_DIR = os.path.join("trees/tree_results/", BASENAME, "ortologs_boot")
    PARALOGS_BOOTSTRAP_OUTPUT_DIR = os.path.join("trees/tree_results/", BASENAME, "paralogs_boot")

    # Make ML Trees using ortological sequences without bootstrap.
    make_trees_ortologs(ORTOLOGS_MSA_PATH, ORTOLOGS_OUTPUT_DIR, CPU_CORES, bootstrap=0, num_processes=NUM_PROCESSES)
    merge_results(ORTOLOGS_OUTPUT_DIR, output_file="all_trees.txt")

    # Make ML Trees using ortological sequences with bootstrap (if greater than zero)
    if BOOTSTRAP > 0:
        make_trees_ortologs(ORTOLOGS_MSA_PATH, ORTOLOGS_BOOTSTRAP_OUTPUT_DIR, CPU_CORES, BOOTSTRAP, NUM_PROCESSES)
        merge_results(ORTOLOGS_BOOTSTRAP_OUTPUT_DIR, output_file="all_trees_bootstrap.txt", eliminate_trees=True, support_threshold=SUPPORT_THRESHOLD)

    # Make ML Trees using PARALOGS sequences without bootstrap.
    make_trees_paralogs(PARALOGS_MSA_RESULTS, PARALOGS_OUTPUT_DIR, CPU_CORES, bootstrap=0, num_processes=NUM_PROCESSES)
    merge_results(PARALOGS_OUTPUT_DIR, output_file="all_trees.txt")

    # Make ML Trees using PARALOGS sequences with bootstrap (if greater than zero)
    if BOOTSTRAP > 0:
        make_trees_paralogs(PARALOGS_MSA_RESULTS, PARALOGS_BOOTSTRAP_OUTPUT_DIR, CPU_CORES, BOOTSTRAP, NUM_PROCESSES)
        merge_results(PARALOGS_BOOTSTRAP_OUTPUT_DIR, output_file="all_trees_bootstrap.txt", eliminate_trees=True, support_threshold=SUPPORT_THRESHOLD)



