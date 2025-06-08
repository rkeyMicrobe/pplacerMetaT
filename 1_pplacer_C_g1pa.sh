#!/bin/bash

#SBATCH --job-name=outFiles/3.1metatreeC.%j  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 96:00:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/3.1metatreeC.%A.out
#SBATCH --array=1-1

# Print time stamp and other info for output files
date;hostname;pwd

echo "Specifying variables from 1_pplacer_B.sh"
GENE=$1
NCORES=$2
CUTOFF=$3
DATASET=$4
GENE_DIR=$5
SAMP_LIST=$6 # Not needed but included
SAMP_DIR=$7  # Not needed but included
subject=$8

echo "Setting main path and subject file variables"
BASE_DIR="data/genes/${DATASET}/${GENE}"
SAMP_FILE="${subject}.6tr.orfs40.fasta.gz"

echo "searching profile against sequence database for each subject file"
echo "cut-off: $CUTOFF"
module purge; ml hmmer
hmmsearch -A "${BASE_DIR}/${GENE}.${subject}.query.sto" -T ${CUTOFF} --incT ${CUTOFF} --cpu ${NCORES} --tblout "${BASE_DIR}/${GENE}.${subject}.hmm_out.tab" "${BASE_DIR}/${GENE}.hmm" "/blue/b.durham/rebeccakey/2_metaTJiwoon/data/${DATASET}/${SAMP_FILE}"

echo "Aligning query hits to reference alignment"
hmmalign -o "${BASE_DIR}/${GENE}.${subject}.aln.sto" --mapali "${BASE_DIR}/${GENE}.sto" "${BASE_DIR}/${GENE}.hmm" "${BASE_DIR}/${GENE}.${subject}.query.sto"

echo "Converting alignment file format from Stockholm back to Fasta using seqmagick"
module purge; ml seqmagick
seqmagick convert "${BASE_DIR}/${GENE}.${subject}.aln.sto" "${BASE_DIR}/${GENE}.${subject}.aln.fasta"

# Run pplacer to place sequences onto a fixed reference tree!
echo "Running pplacer"
module purge; ml pplacer

MAXPEND=0.7
echo "Max-pend value: $MAXPEND"
pplacer -c "${BASE_DIR}/${GENE}.refpkg" --keep-at-most 1 --max-pend $MAXPEND "${BASE_DIR}/${GENE}.${subject}.aln.fasta" -o "results/${DATASET}/${GENE}/${GENE}.${subject}.aln.jplace"

