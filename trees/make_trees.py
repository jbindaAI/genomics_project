from concurrent.futures import ProcessPoolExecutor
import subprocess
import argparse
from tqdm import tqdm
import os


def run_tree_computation(msa_file: str, msa_path: str, output_dir: str, cpu_cores: int, bootstrap: int):
    exact_filepath = os.path.join(msa_path, msa_file)
    output_prefix = os.path.join(output_dir, f"{msa_file.split("-")[0]}")
    
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
            command + ["-b", str(bootstrap)]
        
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


def make_trees(msa_path: str, output_dir: str, cpu_cores: int, bootstrap: int, num_processes: int = 4):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for performing tree computations.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    parser.add_argument("--cpu_cores", required=True, type=int, help="Number of CPU cores to use during computations.")
    parser.add_argument("--num_processes", required=True, type=int, help="Number of separate processes to run. CPU cores are uniformly distributed on processes.")
    parser.add_argument("--bootstrap", required=True, type=int, help="Number of bootstrap replicates.")
    args = parser.parse_args()

    CPU_CORES = args.cpu_cores
    NUM_PROCESSES = args.num_processes
    BOOTSTRAP = args.bootstrap
    BASENAME = args.basename
    MSA_PATH = os.path.join("allignment/msa_results", BASENAME)
    OUTPUT_DIR = os.path.join("trees/tree_results/", BASENAME)

    make_trees(MSA_PATH, OUTPUT_DIR, CPU_CORES, BOOTSTRAP, NUM_PROCESSES)
