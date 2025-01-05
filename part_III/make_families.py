import os
import pickle
import argparse
from collections import defaultdict


def parse_clusters(cluster_file):
    """
    Parse MMseqs2 cluster results to map sequences to clusters.
    
    Parameters:
    - cluster_file: str, Path to cluster_results.tsv file.
    
    Returns:
    - cluster_map: dict, Cluster ID -> List of sequences.
    """
    cluster_map = defaultdict(list)
    with open(cluster_file, "r") as f:
        for line in f:
            cluster, sequence = line.strip().split("\t")
            cluster_map[cluster].append(sequence)
    return cluster_map


def filter_one_to_one_clusters(cluster_map, min_cluster_size, genome_map):
    """
    Filter clusters to include only 1-1 clusters (one sequence per genome).
    
    Parameters:
    - cluster_map: dict, Cluster ID -> List of sequence IDs.
    - genome_map: dict, Sequence ID -> (Genome ID, Genome name, sequence) mapping.
    
    Returns:
    - filtered_clusters: dict, Filtered clusters with genome IDs as keys.
    """
    filtered_clusters = {}
    for cluster, sequences in cluster_map.items():
        if len(sequences) < min_cluster_size:
            continue
        genome_counts = defaultdict(int)
        genome_sequences = {}
        
        for seq in sequences:
            genomeID, _, _ = genome_map.get(seq, "Unknown")
            genome_counts[genomeID] += 1
            genome_sequences[genomeID] = seq
        
        # Check if cluster is 1-1 (exactly one sequence per genome)
        if len(genome_counts) == len(genome_sequences) and all(count == 1 for count in genome_counts.values()):
            filtered_clusters[cluster] = {genomeID: genome_sequences[genomeID] for genomeID in genome_sequences}
    return filtered_clusters


def load_genome_map(dataset:str):
    """
    Loads previously computed genome maps.
    """
    with open(f"part_I/data/maps/{dataset}/genome_map.pkl", "rb") as f1: 
        genome_map = pickle.load(f1)
    return genome_map


def save_filtered_clusters(filtered_clusters, output_file):
    """
    Save filtered 1-1 clusters to a file.
    
    Parameters:
    - filtered_clusters: dict, Filtered clusters with genome names as keys.
    - output_file: str, Path to save the output file.
    """
    with open(output_file, "w") as f:
        for cluster, genomes in filtered_clusters.items():
            cluster_line = f"{cluster}:" + "".join(f"{seq} " for seq in genomes.values())
            f.write(cluster_line + "\n")


def prepare_families(filtered_clusters:dict, genome_map:dict, output_dir:str):
    """
    Prepare protein families for multiple sequence analysis.
    """
    for cluster, genomes in filtered_clusters.items():
        # Save sequences to a file
        with open(os.path.join(output_dir, f"{cluster}.fasta"), "w") as f:
            for genome, seq_ID in genomes.items():
                f.write(f">{genome}\n{genome_map[seq_ID][2]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for filtering 1-1 clusters from MMseqs2 results.")
    parser.add_argument("--accession_filename", required=True, help="Name of an accession file.")
    parser.add_argument("--min_cluster_size", required=True, type=int, help="Minimum number of sequences in clusters.")
    args = parser.parse_args()

    # Paths and parameters
    BASENAME = args.accession_filename.split(".")[0]
    MIN_CLUSTER_SIZE = args.min_cluster_size
    CLUSTER_RES_PATH = os.path.join("part_II/clustering_results", BASENAME, "clustering_results_cluster.tsv")
    FILTERED_OUTPUT_FILE = os.path.join("part_III/filtered_clusters", BASENAME, "filtered_clusters_1-1.txt")
    FAMILIES_OUTPUT_DIR = os.path.join("part_III/protein_families", BASENAME)

    os.makedirs(os.path.join("part_III/filtered_clusters", BASENAME), exist_ok=True)
    os.makedirs(FAMILIES_OUTPUT_DIR, exist_ok=True)

    # Load genome map
    genome_map = load_genome_map(BASENAME)
    print(f"Loaded genome map with {len(genome_map)} sequences.")

    # Parse clusters
    cluster_map = parse_clusters(CLUSTER_RES_PATH)
    print(f"Parsed {len(cluster_map)} clusters from {CLUSTER_RES_PATH}.")

    # Filter 1-1 clusters
    filtered_clusters = filter_one_to_one_clusters(cluster_map, MIN_CLUSTER_SIZE, genome_map)
    print(f"Filtered to {len(filtered_clusters)} 1-1 clusters.")

    # Save filtered clusters in a file
    save_filtered_clusters(filtered_clusters, FILTERED_OUTPUT_FILE)
    print(f"Filtered clusters saved to {FILTERED_OUTPUT_FILE}.")

    # Prepare families for MSA:
    prepare_families(filtered_clusters, genome_map, FAMILIES_OUTPUT_DIR)
    print("Families prepared for multiple sequence analysis.")