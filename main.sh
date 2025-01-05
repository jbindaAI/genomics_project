#!/bin/bash

# Exit script on error
set -e

# Step 0: Set global variables
ACCESSION_FILE="bacteria.txt"

## Clustering options
MIN_SEQ_ID=0.5 # Minimum sequence identity for clustering.
COVERAGE=0.8 # Minimum coverage for clustering.

# Step 1: Prepare data - download proteomes using accessions IDs defined in the ACCESSION_FILE.
echo "Step 1: Downloading proteomes..."
python part_I/prepare_data.py --accession_file "$ACCESSION_FILE" 

# Step 2: Perform clustering with MMseqs2
echo "Step 2: Clustering protein sequences with MMseqs2..."
python part_II/cluster.py --accession_file "$ACCESSION_FILE" --min_seq_id $MIN_SEQ_ID --coverage $COVERAGE

# Step 3: Analyze clusters and extract families (1-to-1)
#echo "Step 4: Analyzing clusters to extract gene families..."
#python analyze_clusters.py --cluster_file "$CLUSTER_DIR/cluster_results.tsv" --output_file "$OUTPUT_PREFIX/family_1to1.tsv"

# Step 4: Multi-sequence alignment

# Step 5: Construct gene trees

# Step 6: Combine gene trees into genome tree

echo "Pipeline completed successfully!"
