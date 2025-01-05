import os
import subprocess
import argparse

def mmseqs2_cluster(combined_fasta_path, output_dir, tmp_dir, min_seq_id=0.5, coverage=0.8):
    """
    Perform clustering of protein sequences using MMseqs2.
    
    Parameters:
    - input_dir: str, Directory containing protein FASTA files.
    - output_dir: str, Path to save MMseqs2 results.
    - tmp_dir: str, Path to temporary directory.
    - min_seq_id: float, Minimum sequence identity for clustering.
    - coverage: float, Minimum coverage for clustering.
    """

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Change working directory to the desired output directory. (MMseqs2 will save results there)
    original_dir = os.getcwd()
    os.chdir(output_dir)

    combined_fasta_path = os.path.join(original_dir, combined_fasta_path)
    tmp_dir = os.path.join(original_dir, tmp_dir)

    try:
        output_name = "clustering_results"
        # Run MMseqs2 easy-cluster
        command = [
            "mmseqs", "easy-cluster",
            combined_fasta_path, output_name, tmp_dir,
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage)
        ]

        print(f"Running MMseqs2 clustering: {' '.join(command)}")
        subprocess.run(command, check=True)
        print(f"Clustering complete.")
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for performing clustering with MMSeqs2 algorithm.")
    parser.add_argument("--accession_filename", required=True, help="Name of an accession file.")
    parser.add_argument("--min_seq_id", required=True, type=float, help="Minimum sequence identity for clustering.")
    parser.add_argument("--coverage", required=True, type=float, help="Minimum coverage for clustering.")
    args = parser.parse_args()

    BASENAME = args.accession_filename.split(".")[0]
    COMBINED_FASTA_PATH = os.path.join("part_I/data/combined_fasta", BASENAME, "combined_proteins.faa")
    OUTPUT_DIR = os.path.join("part_II/clustering_results/", BASENAME)        # Output name for MMseqs2 results
    TMP_DIR = "part_II/tmp"                   # Temporary directory for MMseqs2
    
    os.makedirs(TMP_DIR, exist_ok=True)

    mmseqs2_cluster(COMBINED_FASTA_PATH, OUTPUT_DIR, TMP_DIR, min_seq_id=args.min_seq_id, coverage=args.coverage)
