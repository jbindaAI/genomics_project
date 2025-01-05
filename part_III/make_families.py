import os
import pickle
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
            sequence, cluster = line.strip().split("\t")
            cluster_map[cluster].append(sequence)
    return cluster_map


def filter_one_to_one_clusters(cluster_map, genome_map):
    """
    Filter clusters to include only 1-1 clusters (one sequence per genome).
    
    Parameters:
    - cluster_map: dict, Cluster ID -> List of sequence IDs.
    - genome_map: dict, Sequence -> Genome mapping from sequence ID to genome ID.
    
    Returns:
    - filtered_clusters: dict, Filtered clusters with genome IDs as keys.
    """
    filtered_clusters = {}
    for cluster, sequences in cluster_map.items():
        if len(sequences) < 2:
            continue
        genome_counts = defaultdict(int)
        genome_sequences = {}
        
        for seq in sequences:
            genome = genome_map.get(seq, "Unknown")
            genome_counts[genome] += 1
            genome_sequences[genome] = seq
        
        # Check if cluster is 1-1 (exactly one sequence per genome)
        if len(genome_counts) == len(genome_sequences) and all(count == 1 for count in genome_counts.values()):
            filtered_clusters[cluster] = {genome: genome_sequences[genome] for genome in genome_sequences}
    
    return filtered_clusters


def load_genome_maps(dataset:str):
    """
    Loads previously computed genome maps.
    """
    with open(f"part_I/data/maps/{dataset}_seqID2genomeID.pkl", "rb") as f1:
        seqID2genomeID = pickle.load(f1)

    with open(f"part_I/data/maps/{dataset}_genomeID2name.pkl", "rb") as f2: 
        genomeID2name = pickle.load(f2)

    return seqID2genomeID, genomeID2name


def save_filtered_clusters(filtered_clusters, output_file):
    """
    Save filtered 1-1 clusters to a file.
    
    Parameters:
    - filtered_clusters: dict, Filtered clusters with genome names as keys.
    - output_file: str, Path to save the output file.
    """
    with open(output_file, "w") as f:
        for cluster, genomes in filtered_clusters.items():
            cluster_line = f"{cluster}:\t" + "\t".join(seq for seq in genomes.values())
            f.write(cluster_line + "\n")


if __name__ == "__main__":
    # Paths and parameters
    cluster_file = "part_II/ClusterRes_cluster.tsv"
    protein_files_dir = "part_I/data/proteome_fasta/"
    output_file = "part_III/filtered_clusters_1-1.txt"
    dataset = "bacteria"

    # Load genome maps
    seqID2genomeID, genomeID2name = load_genome_maps(dataset)
    print(f"Loaded genome map with {len(seqID2genomeID)} sequences.")

    # Parse clusters
    cluster_map = parse_clusters(cluster_file)
    print(f"Parsed {len(cluster_map)} clusters from {cluster_file}.")

    # Filter 1-1 clusters
    filtered_clusters = filter_one_to_one_clusters(cluster_map, seqID2genomeID)
    print(f"Filtered to {len(filtered_clusters)} 1-1 clusters.")

    # Save filtered clusters
    save_filtered_clusters(filtered_clusters, output_file)
    print(f"Filtered clusters saved to {output_file}.")