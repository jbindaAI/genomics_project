import subprocess
import argparse
import os


def make_consensus_tree(all_trees_path: str, output_dir: str, min_support: float, cpu_cores: int):
    os.makedirs(output_dir, exist_ok=True)

    output_prefix = os.path.join(output_dir, "consensus_tree")

    try:
        # Prepare the IQ-TREE command
        command = [
            "iqtree",
            "-con", all_trees_path,
            "-minsup", str(min_support),
            "-nt", str(cpu_cores),
            "-pre", output_prefix,
            #"-quiet"
        ]
        
        # Run the IQ-TREE command
        result = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Return success
        return ("Success")
    except subprocess.CalledProcessError as e:
        return (f"Failed (Code {e.returncode})")
    except FileNotFoundError as e:
        return ("File Not Found")
    except Exception as e:
        return (f"Unexpected Error: {e}")
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script computes a consensus tree.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    parser.add_argument("--min_support", required=True, type=float, help="Value of Consensus support threshold.")
    parser.add_argument("--cpu_cores", required=True, type=int, help="NUmber of CPU cores.")
    args = parser.parse_args()

    BASENAME = args.basename
    MIN_SUPPORT = args.min_support
    CPU_CORES = args.cpu_cores
    ALL_TREES_PATH = os.path.join("trees/tree_results/", BASENAME, "unrooted", "all_trees.txt")
    ALL_TREES_PATH_BOOTSTRAP = os.path.join("trees/tree_results/", BASENAME, "unrooted_boot", "all_trees_bootstrap.txt")

    OUTPUT_DIR = os.path.join("trees/consensus_results", BASENAME)

    # Make Consensus Tree from ML Trees build on Ortological sequences WITHOUT Bootstrap
    make_consensus_tree(ALL_TREES_PATH, os.path.join(OUTPUT_DIR, "vanilla"), MIN_SUPPORT, CPU_CORES)

    # Make Consensus Tree from ML Trees build on Ortological sequences WITH Bootstrap
    make_consensus_tree(ALL_TREES_PATH_BOOTSTRAP, os.path.join(OUTPUT_DIR, "boostrapped"), MIN_SUPPORT, CPU_CORES)


