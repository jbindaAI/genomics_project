#!/bin/bash

# Exit script on error
set -e

# Step 0: Set global variables
ACCESSION_FILE="bacteria.txt" # PUT THERE NAME OF ACCESION FILE

## Clustering options
MIN_SEQ_ID=0.5      # Minimum sequence identity for clustering.
COVERAGE=0.8        # Minimum coverage for clustering.
MIN_CLUSTER_SIZE=4  # Minimum number of sequences in clusters.

BASENAME="${ACCESSION_FILE%.*}"  # Removes the extension, used for saving intermediate results.

# Step 1: Prepare data - download proteomes using accessions IDs defined in the ACCESSION_FILE.
echo "Step 1: Downloading proteomes..."
python data_preparation/prepare_data.py --accession_file "$ACCESSION_FILE" 

# Step 2: Perform clustering with MMseqs2
echo "Step 2: Clustering protein sequences with MMseqs2..."
python clustering/cluster.py --basename "$BASENAME" --min_seq_id $MIN_SEQ_ID --coverage $COVERAGE

# Step 3: Analyze clusters and extract families (1-to-1)
echo "Step 4: Analyzing clusters to extract gene families..."
python families/make_families.py --basename "$BASENAME" --min_cluster_size $MIN_CLUSTER_SIZE

# Step 4: Multi-sequence alignment
echo "Step 5: Performing multiple sequence alignments..."
python allignment/allign.py --basename "$BASENAME"

# Step 5: Construct gene trees

# Step 6: Combine gene trees into genome tree

echo "Pipeline completed successfully!"
