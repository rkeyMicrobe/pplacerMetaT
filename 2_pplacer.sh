#!/bin/bash

#SBATCH --job-name=outFiles/4.2metatree.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 0:10:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/4.2metatree.%A.out
#SBATCH --array=1-26

# Print out the date, hostname and working directory
date;hostname;pwd

echo "Setting variables"
GENE=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genesAbbr.txt)
GENE_DIR=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genes.txt)
DATASET="g1pa"

echo "Generating a base tree with all jplace information"
module purge; ml pplacer
guppy fat -o "results/${DATASET}/${GENE}/${GENE}.${DATASET}.all.xml" --average results/${DATASET}/${GENE}/${GENE}.*.jplace

echo "Running the treeColor_eukNprok.py script to place taxa colors and count information (edge widths) to base tree"
module purge; ml python/2.7
SCRIPT="tools/treeColor_eukNprok/treeColor_eukNprok.py"
python ${SCRIPT} -o "results/${DATASET}/${GENE}/${GENE}.${DATASET}.allfatcolor.xml" "results/${DATASET}/${GENE}/${GENE}.${DATASET}.all.xml"
python ${SCRIPT} -o "results/${DATASET}/${GENE}/${GENE}.${DATASET}.allcolor.xml" -a "results/${DATASET}/${GENE}/${GENE}.${DATASET}.all.xml"
    
echo "Generating list of sample list from jplaces for loop below"
ls results/${DATASET}/${GENE}/${GENE}.*.jplace | sort | uniq > "results/${DATASET}/${GENE}/${GENE}.${DATASET}.sample_list.txt"

echo "Processing CSVs by summed samples"
module purge; ml pplacer   
for sample in $(cat "results/${DATASET}/${GENE}/${GENE}.${DATASET}.sample_list.txt"); do
    guppy to_csv -o "${sample}.taxID.csv" "${sample}"
done

echo "Moving taxID files to their own directory"
mv results/${DATASET}/${GENE}/*.taxID.csv results/${DATASET}/${GENE}/${GENE}_taxIDs

echo "Compress all .jplace files into one tarball zip"
tar czf results/${DATASET}/${GENE}/${GENE}.jplace.tar.gz results/${DATASET}/${GENE}/${GENE}.*.jplace
rm results/${DATASET}/${GENE}/${GENE}.*.jplace
