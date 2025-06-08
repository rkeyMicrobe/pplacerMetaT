#!/bin/bash

#SBATCH --job-name=outFiles/5.G2NSmORF.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 1:00:00
#SBATCH --output=outFiles/5.G2NSmORF.%A.out
#SBATCH --array=1-26

# Print date, hostname and working directory
date;hostname;pwd

# Define variables
echo "Defining variables"
GENE=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genesAbbr.txt)
GENE_DIR=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genes.txt)
DATASET="g1pa"
INGROUP="data/genes/edges/${GENE}.just_edges.txt"
TREE_COL="tools/treeColor_eukNprok/treeColor_ids/treeColor_key.csv"

##################################################################################3

echo "Beginning final step of analysis..."

NORMFACTORS="data/normalize/${DATASET}_normFactors.csv"
SCRIPT="tools/treeColor_eukNprok/${DATASET}_pplacer_idClassify_counts.py" 

echo "Norm factor file: $NORMFACTORS"
echo "script file: $SCRIPT"

echo "Make list of csv files for loop"
CSV_FILE=$(ls results/${DATASET}/${GENE}/${GENE}_taxIDs/${GENE}.*.taxID.csv | head -1)

echo "Running pplacer_idClassify_counts script"
module purge; ml python/2.7
python ${SCRIPT} -e -g -i ${INGROUP} -c ${TREE_COL} -n ${NORMFACTORS}  ${CSV_FILE} > "results/${DATASET}/${GENE}/${GENE}.${DATASET}.counts.csv"

# Process each csv file
echo "Processing each taxID.csv"
for CSV_FILE in $(ls results/${DATASET}/${GENE}/${GENE}_taxIDs/${GENE}.*.taxID.csv); do 
    echo "Processing ${CSV_FILE}..."
    python ${SCRIPT} -g -i ${INGROUP} -c ${TREE_COL} -n ${NORMFACTORS} ${CSV_FILE} >> "results/${DATASET}/${GENE}/${GENE}.${DATASET}.counts.csv"; 
done 
