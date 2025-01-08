#!/usr/bin/env bash
#DIR="./genomika"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate maft

for file in *fasta; do
        echo $(pwd)
        mafft --auto "$file" > "${file}_msa"
done
