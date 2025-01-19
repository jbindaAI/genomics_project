from typing import Literal
import subprocess
import argparse
import os

def run_r_script(tree_file, output_file, method:Literal["MRP", "RF", "SPR"], n_cores:int):
    rscript_path = "Rscript" #Needed to run scripts written in R
    r_script = "trees/super_tree.R"
    multicore = "True" if n_cores>1 else "False"
    try:
        subprocess.run([rscript_path, r_script, tree_file, output_file, method, multicore, str(n_cores)], check=True)
        print(f"Supertree created successfully: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running R script: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script computes a consensus tree.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    parser.add_argument("--method", required=True, help="Method used to infer SuperTree. Available: MRP, RF, SPR")
    parser.add_argument("--cpu_cores", required=True, type=int, help="Number of CPU cores.")
    args = parser.parse_args()

    BASENAME = args.basename
    METHOD = args.method
    CPU_CORES = args.cpu_cores
    ALL_TREES_PATH = os.path.join("trees/tree_results/", BASENAME, "ortologs", "all_trees.txt")
    ALL_TREES_PATH_BOOTSTRAP = os.path.join("trees/tree_results/", BASENAME, "ortologs_boot", "all_trees_bootstrap.txt")

    OUTPUT_DIR = os.path.join("trees/super_tree_results/ortologs", BASENAME)

    # Make SuperTree from ML Trees build on Paralogical sequences WITHOUT Bootstrap
    os.makedirs(os.path.join(OUTPUT_DIR, "vanilla"), exist_ok=True)
    run_r_script(ALL_TREES_PATH, os.path.join(OUTPUT_DIR, "vanilla", "supertree.nwk"), METHOD, CPU_CORES)

    # Make SuperTree from ML Trees build on Paralogical sequences WITH Bootstrap
    os.makedirs(os.path.join(OUTPUT_DIR, "boostrapped"), exist_ok=True)
    run_r_script(ALL_TREES_PATH_BOOTSTRAP, os.path.join(OUTPUT_DIR, "boostrapped", "supertree.nwk"), METHOD, CPU_CORES)

