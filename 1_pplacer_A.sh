#!/bin/bash

#SBATCH --job-name=outFiles/1.1metatreeA.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mem=2GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 00:15:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/1.%a.1metatreeA.%A.out
#SBATCH --array=1-16

# Print time stamp and other info for output files
date;hostname;pwd

echo "Set variables"
# Variables
GENE=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genesAbbr.txt)
GENE_DIR=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/genes/genes.txt)
CUTOFF="40"
NCORES="4"
DATASET="g2pa"
# Variable paths
TAX_DB="data/taxonomy.db"
BASE_DIR="data/genes/genes_start/${GENE_DIR}"
GENE_DIR="data/genes/${DATASET}/${GENE_DIR}"

# Python
ml python
echo "Start python"
python tools/make_seq_info_all.py ${BASE_DIR}/${GENE}.trimal2.fasta -o ${GENE_DIR}/${GENE}.seq_info_all.csv

# Taxtastic
module purge; ml taxtastic
echo "Start taxit"
taxit update_taxids -o ${GENE_DIR}/${GENE}.seq_info_all.updated.csv ${GENE_DIR}/${GENE}.seq_info_all.csv $TAX_DB 

# Create taxtable
awk -F, '{print $2}' ${GENE_DIR}/${GENE}.seq_info_all.updated.csv | sort | uniq | sed '/tax_id/d' | tr -d '"' > ${GENE_DIR}/${GENE}.tax_ids_all.txt
taxit taxtable $TAX_DB -f ${GENE_DIR}/${GENE}.tax_ids_all.txt -o ${GENE_DIR}/${GENE}.taxa.csv

# Make refpkg
echo "Making the refpkg on $GENE in $GENE_DIR"
taxit create -l $GENE -P ${GENE_DIR}/${GENE}.refpkg \
--taxonomy ${GENE_DIR}/${GENE}.taxa.csv \
--aln-fasta ${BASE_DIR}/${GENE}.trimal2.fasta \
--seq-info ${GENE_DIR}/${GENE}.seq_info_all.updated.csv \
--tree-stats ${BASE_DIR}/RAxML_info.brL.${GENE}.tre \
--tree-file ${BASE_DIR}/RAxML_result.brL.${GENE}.tre \
--no-reroot

# Seqmagick
module purge; ml seqmagick
echo "convert fasta format to Stockholm format for hmmer package use"
seqmagick convert ${BASE_DIR}/${GENE}.trimal2.fasta ${GENE_DIR}/${GENE}.sto

# Hmmer
module purge; ml hmmer
echo "build profile of multiple alignment from seq info stored in .sto file"
hmmbuild ${GENE_DIR}/${GENE}.hmm ${GENE_DIR}/${GENE}.sto

# Move to Script B
echo "Move each gene hmm alignment from the array to the Script B"
 sbatch 1_pplacer_B.sh $GENE $NCORES $CUTOFF $DATASET $GENE_DIR

