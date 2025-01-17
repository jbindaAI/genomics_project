import subprocess
import argparse
from tqdm import tqdm
import os
from concurrent.futures import ProcessPoolExecutor


def run_mafft(fasta_file: str, fasta_path: str, output_dir: str):
    """
    Run MAFFT alignment on a single fasta file.
    """
    exact_filepath = os.path.join(fasta_path, fasta_file)
    output_file = os.path.join(output_dir, f"{os.path.splitext(fasta_file)[0]}-aligned.fasta")
    try:
        command = [
            "mafft",
            "--quiet",
            "--auto",
            exact_filepath,
        ]
        with open(output_file, "w") as output:
            subprocess.run(command, check=True, stdout=output, text=True)
        return (fasta_file, "Success")
    except subprocess.CalledProcessError as e:
        return (fasta_file, f"Failed (Code {e.returncode})")
    except FileNotFoundError as e:
        return (fasta_file, "File Not Found")
    except Exception as e:
        return (fasta_file, f"Unexpected Error: {e}")


def mafft_align(fasta_path: str, output_dir: str, num_processes: int = 4):
    """
    Perform MAFFT alignment in parallel.
    """
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = os.listdir(fasta_path)
    
    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # Wrap the executor map with tqdm for progress tracking
        results = list(tqdm(
            executor.map(
                run_mafft,
                fasta_files,
                [fasta_path] * len(fasta_files),
                [output_dir] * len(fasta_files),
            ),
            total=len(fasta_files),
            desc="Aligning sequences",
            unit="file"
        ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for performing multiple sequence alignment with MAFFT algorithm.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    parser.add_argument("--num_processes", type=int, default=4, help="Number of processes to use for multiprocessing.")
    args = parser.parse_args()

    BASENAME = args.basename
    NUM_PROCESSES = args.num_processes

    PROTEIN_FAMILIES_PATH_ORTOLOGS = os.path.join("families/protein_families/ortologs", BASENAME)
    PROTEIN_FAMILIES_PATH_PARALOGS = os.path.join("families/protein_families/paralogs", BASENAME)

    OUTPUT_DIR_ORTOLOGS = os.path.join("allignment/msa_results/ortologs", BASENAME)
    OUTPUT_DIR_PARALOGS = os.path.join("allignment/msa_results/paralogs", BASENAME)

    # Perform MAFFT alignment in parallel
    print("Preparing MSA for ortological sequences ...")
    mafft_align(PROTEIN_FAMILIES_PATH_ORTOLOGS, OUTPUT_DIR_ORTOLOGS, NUM_PROCESSES)
    print("Done.")
    print("Preparing MSA for paralogical sequences ...")
    mafft_align(PROTEIN_FAMILIES_PATH_PARALOGS, OUTPUT_DIR_PARALOGS, NUM_PROCESSES)
    print("Done.")
