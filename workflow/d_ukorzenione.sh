DIR=" "
INPUTS="${DIR}/msa_fasta"
OUTPUTS="${DIR}/fast_trees"

for file in ${INPUTS}/*fasta
do
        echo $file
        NAME=$(basename "$file" .fasta)
        ./FastTree < "$file" > "${OUTPUTS}/${NAME}.nwk"
done
