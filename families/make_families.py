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


def filter_clusters(cluster_map, min_cluster_size, genome_map, one2one:bool):
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
            genomeID, _, _ = genome_map.get(seq, "Unknown")
            genome_counts[genomeID] += 1
            genome_sequences[genomeID] = seq
        
        # Check if cluster is 1-1 (exactly one sequence per genome and we want only clusters with ALL sequences - bijective)
        if one2one and len(genome_counts) == len(genome_sequences) and len(genome_counts)==NUMBER_OF_ALL_SEQUENCES and all(count == 1 for count in genome_counts.values()):
            filtered_clusters[cluster] = {genomeID: genome_sequences[genomeID] for genomeID in genome_sequences}
        elif not one2one:
            filtered_clusters[cluster] = {genomeID: genome_sequences[genomeID] for genomeID in genome_sequences}

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

    os.makedirs(CLUSTER_OUTPUT_DIR_ORTOLOGS, exist_ok=True)
    os.makedirs(CLUSTER_OUTPUT_DIR_PARALOGS, exist_ok=True)

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
    clusters_paralogs = filter_clusters(cluster_map, MIN_CLUSTER_SIZE, genome_map, one2one=False)
    print(f"Extracted {len(clusters_paralogs)} clusters with paralogs.")

    # Filter 1-1 clusters
    clusters_ortologs = filter_clusters(cluster_map, MIN_CLUSTER_SIZE, genome_map, one2one=True)
    print(f"Extracted {len(clusters_ortologs)} 1-1 clusters (without paralogs, bijective).")

    # Save clusters with paralogs in a file
    paralogs_save_path = os.path.join(CLUSTER_OUTPUT_DIR_PARALOGS, "clusters_paralogs.txt")
    save_filtered_clusters(clusters_paralogs, paralogs_save_path)
    print(f"Clusters with paralogs saved to {paralogs_save_path}.")

    # Save clusters without paralogs in a file
    ortologs_save_path = os.path.join(CLUSTER_OUTPUT_DIR_ORTOLOGS, "clusters_1-1.txt")
    save_filtered_clusters(clusters_ortologs, ortologs_save_path)
    print(f"Clusters without paralogs saved to {ortologs_save_path}.")

    # Prepare families for MSA (with paralogs):
    prepare_families(clusters_paralogs, genome_map, PARALOGS_FAMILIES_OUTPUT_DIR)

    # Prepare families for MSA (without paralogs):
    prepare_families(clusters_ortologs, genome_map, ORTOLOGS_FAMILIES_OUTPUT_DIR)
    print("Families prepared for multiple sequence analysis.")