#!/bin/bash

# Exit script on error
set -e

# Step 0: Set global variables
ACCESSION_FILE="bacteria.txt" # PUT THERE NAME OF ACCESION FILE

## Clustering options
MIN_SEQ_ID=0.5      # Minimum sequence identity for clustering.
COVERAGE=0.8        # Minimum coverage for clustering.
MIN_CLUSTER_SIZE=4  # Minimum number of sequences in clusters.

## MSA options
MSA_NUM_PROCESSES=4

## Tree options
CPU_CORES=12                      # Number of CPU cores to use during tree computation (If you don't know, try: os.cpu_count())
TREE_NUM_PROCESSES=4              # Number of separate processes to run. Note, it would be better if: CPU_CORES % NUM_PROCESSES == 0
BOOTSTRAP_REPLICATES=10           # Number of bootstrap replicates. If greater than zero, trees will be computed two times. Once without bootstraping and second one with apllying bootstrap.
BOOTSTRAP_SUPPORT_THRESHOLD=70.0  # Bootstrap trees with mean support lower than threshold are eliminated from analysis.

## Consensus Tree options
MIN_SUPPORT=0           # Value from 0 to 1. If zero it perform Greedy Consensus, if 0.5 it performs Majority Consensus
CONSENSUS_CPU_CORES=4


BASENAME="${ACCESSION_FILE%.*}"

# Step 1: Prepare data - download proteomes using accessions IDs defined in the ACCESSION_FILE.
echo "Step 1: Downloading proteomes..."
python3 data_preparation/prepare_data.py --accession_file "$ACCESSION_FILE" 

# # Step 2: Perform clustering with MMseqs2
echo "Step 2: Clustering protein sequences with MMseqs2..."
python3 clustering/cluster.py --basename "$BASENAME" --min_seq_id $MIN_SEQ_ID --coverage $COVERAGE

# # Step 3: Analyze clusters and extract families (1-to-1)
echo "Step 4: Analyzing clusters to extract gene families..."
python3 families/make_families.py --basename "$BASENAME" --min_cluster_size $MIN_CLUSTER_SIZE

# # Step 4: Multi-sequence alignment
echo "Step 4: Performing multiple sequence alignments..."
python3 allignment/allign.py --basename "$BASENAME" --num_processes "$MSA_NUM_PROCESSES"

# # Step 5: Construct gene trees
echo "Step 5: Constructing family trees..."
python3 trees/make_trees.py --basename "$BASENAME" --cpu_cores "$CPU_CORES" --bootstrap "$BOOTSTRAP_REPLICATES" --num_processes "$TREE_NUM_PROCESSES" --support_threshold "$BOOTSTRAP_SUPPORT_THRESHOLD"

# Step 6: Construct Consensus Trees
echo "Step 6: Constructing Consensus trees..."
python3 trees/make_consensus_tree.py --basename "$BASENAME" --min_support "$MIN_SUPPORT" --cpu_cores "$CONSENSUS_CPU_CORES"

echo "Pipeline completed successfully!"
