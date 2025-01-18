import os
import pickle
import argparse
from collections import defaultdict


def parse_clusters(cluster_file: str):
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


def filter_clusters_ortologs(cluster_map, min_cluster_size, genome_map):
    """
    Filter clusters to include cluster of size >= min_cluster_size.
    Also If needed select only 1-1 clusters (one sequence per genome).
    
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
            _, genome_name, _ = genome_map.get(seq, "Unknown") # SequenceID : (Genome ID, Genome name, sequence)
            genome_counts[genome_name] += 1
            genome_sequences[genome_name] = seq
        
        # Check if cluster is 1-1 (exactly one sequence per genome and we want only clusters with ALL sequences - bijective)
        if len(genome_counts) == len(genome_sequences) and len(genome_counts)==NUMBER_OF_ALL_SEQUENCES and all(count == 1 for count in genome_counts.values()):
            filtered_clusters[cluster] = {genome_name: genome_sequences[genome_name] for genome_name in genome_sequences}

    return filtered_clusters


def filter_clusters_paralogs(cluster_map, min_cluster_size):
    """Filter out clusters of size lower than threshold"""
    filtered_clusters = {}
    for cluster, sequences in cluster_map.items():
        if len(sequences) < min_cluster_size:
            continue

        filtered_clusters[cluster] = sequences

    return filtered_clusters


def load_genome_map(dataset:str):
    """
    Loads previously computed genome maps.
    """
    with open(f"data_preparation/data/maps/{dataset}/genome_map.pkl", "rb") as f1: 
        genome_map = pickle.load(f1)
    
    with open(f"data_preparation/data/maps/{dataset}/genomeID2name.pkl", "rb") as f2:
        genomeID2name = pickle.load(f2)
    
    return genome_map, genomeID2name


def prepare_ortologs_families(filtered_clusters:dict, genome_map:dict, output_dir:str):
    """
    Prepare protein families for multiple sequence analysis.
    """
    for cluster, genomes in filtered_clusters.items():
        # Save sequences to a file
        with open(os.path.join(output_dir, f"{cluster}.fasta"), "w") as f:
            for genome, seq_ID in genomes.items():
                f.write(f">{genome}\n{genome_map[seq_ID][2]}\n")


def prepare_paralogs_families(filtered_clusters:dict, genome_map:dict, output_dir:str):
    """
    Prepare protein families for multiple sequence analysis.
    """
    for cluster, seq_ids in filtered_clusters.items():
        # Save sequences to a file
        with open(os.path.join(output_dir, f"{cluster}.fasta"), "w") as f:
            for seq_ID in seq_ids:
                genome_name = genome_map[seq_ID][1]
                f.write(f">{genome_name}\n{genome_map[seq_ID][2]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for filtering 1-1 clusters from MMseqs2 results.")
    parser.add_argument("--basename", required=True, help="Name used for intermediate results saving.")
    parser.add_argument("--min_cluster_size", required=True, type=int, help="Minimum number of sequences in clusters.")
    args = parser.parse_args()

    # Paths and parameters
    BASENAME = args.basename
    MIN_CLUSTER_SIZE = args.min_cluster_size
    CLUSTER_RES_PATH = os.path.join("clustering/clustering_results", BASENAME, "clustering_results_cluster.tsv")
    
    CLUSTER_OUTPUT_DIR_ORTOLOGS = os.path.join("families/clusters_ortologs", BASENAME)
    CLUSTER_OUTPUT_DIR_PARALOGS = os.path.join("families/clusters_paralogs", BASENAME)

    ORTOLOGS_FAMILIES_OUTPUT_DIR = os.path.join("families/protein_families/ortologs", BASENAME)
    PARALOGS_FAMILIES_OUTPUT_DIR = os.path.join("families/protein_families/paralogs", BASENAME)


    os.makedirs(ORTOLOGS_FAMILIES_OUTPUT_DIR, exist_ok=True)
    os.makedirs(PARALOGS_FAMILIES_OUTPUT_DIR, exist_ok=True)

    # Load genome map
    genome_map, genomeID2name = load_genome_map(BASENAME)
    print(f"Loaded genome map with {len(genome_map)} sequences.")

    NUMBER_OF_ALL_SEQUENCES = len(genomeID2name.keys())

    # Parse clusters
    cluster_map = parse_clusters(CLUSTER_RES_PATH)
    print(f"Parsed {len(cluster_map)} clusters from {CLUSTER_RES_PATH}.")

    # Prepare clusters with paralogs
    clusters_paralogs = filter_clusters_paralogs(cluster_map, MIN_CLUSTER_SIZE)
    print(f"Extracted {len(clusters_paralogs)} clusters with paralogs.")

    # Filter 1-1 clusters
    clusters_ortologs = filter_clusters_ortologs(cluster_map, MIN_CLUSTER_SIZE, genome_map)
    print(f"Extracted {len(clusters_ortologs)} 1-1 clusters (without paralogs, bijective).")

    # Prepare families for MSA (with paralogs):
    prepare_paralogs_families(clusters_paralogs, genome_map, PARALOGS_FAMILIES_OUTPUT_DIR)

    # Prepare families for MSA (without paralogs):
    prepare_ortologs_families(clusters_ortologs, genome_map, ORTOLOGS_FAMILIES_OUTPUT_DIR)