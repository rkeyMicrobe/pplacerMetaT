#!/bin/bash

#SBATCH --job-name=outFiles/make_geneDirectories.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=1GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 00:10:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/1.%a.make_geneDirectories.%A.out
#SBATCH --array=1-26

# Print time stamp and other info for output files
date;hostname;pwd

echo "Set variables"
# Variables
GENE=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genesAbbr.txt)
GENE_DIR=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genes.txt)

DATA_PATH="data/genes/"

echo "making the gene folders w/in data directory"
mkdir ${DATA_PATH}/g1ns/${GENE}
mkdir ${DATA_PATH}/g1pa/${GENE}
mkdir ${DATA_PATH}/g2ns/${GENE}
mkdir ${DATA_PATH}/g2pa/${GENE}


echo "making the result folders and subfolders w/in result directory"
mkdir results/g1ns/${GENE}
mkdir results/g1ns/${GENE}/${GENE}_finals
mkdir results/g1ns/${GENE}/${GENE}_taxIDs

mkdir results/g1pa/${GENE}
mkdir results/g1pa/${GENE}/${GENE}_finals
mkdir results/g1pa/${GENE}/${GENE}_taxIDs

mkdir results/g2ns/${GENE}
mkdir results/g2ns/${GENE}/${GENE}_finals
mkdir results/g2ns/${GENE}/${GENE}_taxIDs

mkdir results/g2pa/${GENE}
mkdir results/g2pa/${GENE}/${GENE}_finals
mkdir results/g2pa/${GENE}/${GENE}_taxIDs


