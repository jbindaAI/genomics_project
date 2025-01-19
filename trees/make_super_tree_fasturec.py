import subprocess
import argparse
import os


def run_fasturec(tree_file: str, output_file: str):
    """
    Removes genome IDs from the tree file using sed and runs Fasturec to generate a supertree.
    """
    try:
        if not os.path.exists(tree_file):
            raise FileNotFoundError(f"Tree file not found: {tree_file}")

        # Use sed to remove genome IDs in the format _[number]
        sed_command = ["sed", "-i", 's/_[0-9]\\+//g', tree_file]
        subprocess.run(sed_command, check=True)

        # Run Fasturec to compute the supertree
        fasturec_command = ["fasturec/fasturec", "-G", tree_file, "-Z", "-e", "-a"]
        subprocess.run(fasturec_command, check=True)

        with open(output_file, "w") as f_out:
            with open('fu.txt') as f_in:
                cost_tree = f_in.readline()
                tree = cost_tree.split(" ")[-1].strip()
                f_out.write(tree)
        
        subprocess.run(["rm", "fu.txt"])

        print(f"Supertree created successfully: {output_file}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script computes a SuperTree.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    args = parser.parse_args()

    BASENAME = args.basename
    ALL_TREES_PATH = os.path.join("trees/tree_results/", BASENAME, "paralogs", "all_trees.txt")
    ALL_TREES_PATH_BOOTSTRAP = os.path.join("trees/tree_results/", BASENAME, "paralogs_boot", "all_trees_bootstrap.txt")

    OUTPUT_DIR = os.path.join("trees/super_tree_results/paralogs", BASENAME)

    # Make SuperTree from ML Trees build on Paralogical sequences WITHOUT Bootstrap
    os.makedirs(os.path.join(OUTPUT_DIR, "vanilla"), exist_ok=True)
    run_fasturec(ALL_TREES_PATH, os.path.join(OUTPUT_DIR, "vanilla", "supertree.nwk"))

    # Make SuperTree from ML Trees build on Paralogical sequences WITH Bootstrap
    os.makedirs(os.path.join(OUTPUT_DIR, "boostrapped"), exist_ok=True)
    run_fasturec(ALL_TREES_PATH_BOOTSTRAP, os.path.join(OUTPUT_DIR, "boostrapped", "supertree.nwk"))

