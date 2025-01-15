import subprocess
import argparse
from tqdm import tqdm
import os


def mafft_align(fasta_path: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_files = os.listdir(fasta_path)
    
    # Initialize tqdm progress bar
    with tqdm(total=len(fasta_files), desc="Aligning sequences", unit="file") as pbar:
        for fasta_file in fasta_files:
            exact_filepath = os.path.join(fasta_path, fasta_file)
            output_file = os.path.join(output_dir, f"{os.path.splitext(fasta_file)[0]}_aligned.fasta")
            try:
                command = [
                    "mafft",
                    "--quiet",
                    "--auto",
                    exact_filepath,
                ]
                with open(output_file, "w") as output:
                    subprocess.run(command, check=True, stdout=output, text=True)
                # Update progress bar on success
                pbar.set_postfix({"Status": "Success"})
            except subprocess.CalledProcessError as e:
                pbar.set_postfix({"Status": f"Failed (Code {e.returncode})"})
                print(f"Error: MAFFT failed for {fasta_file} with return code {e.returncode}. Command: {' '.join(e.cmd)}")
            except FileNotFoundError as e:
                pbar.set_postfix({"Status": "File Not Found"})
                print(f"Error: {e}. Ensure 'mafft' is installed and accessible in your PATH.")
            except Exception as e:
                pbar.set_postfix({"Status": "Unexpected Error"})
                print(f"An unexpected error occurred while processing {fasta_file}: {e}")
            finally:
                # Increment the progress bar after each file
                pbar.update(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for performing multiple sequence allignment with MAFFT algorithm.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    args = parser.parse_args()

    BASENAME = args.basename
    PROTEIN_FAMILIES_PATH = os.path.join("families/protein_families", BASENAME)
    OUTPUT_DIR = os.path.join("allignment/msa_results/", BASENAME)        # Output name for MAFFT results


    mafft_align(PROTEIN_FAMILIES_PATH, OUTPUT_DIR)
