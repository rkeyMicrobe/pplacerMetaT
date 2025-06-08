#!/bin/bash

#SBATCH --job-name=outFiles/2.1metatreeB.%j   # Job name. %j is replaced by the job ID
#SBATCH --ntasks=1                            # The number of tasks to run
#SBATCH --cpus-per-task=1                     # The number of CPUs per task
#SBATCH --nodes=1                             # The number of nodes to use
#SBATCH --mem=2GB                             # The amount of memory to allocate
#SBATCH --account=b.durham                    # The account to use for resources
#SBATCH --qos=b.durham-b                      # Quality of service, specifies priority
#SBATCH -t 00:10:00                           # Time limit in hh:mm:ss
#SBATCH --partition=hpg-milan                 # The partition (queue) to submit to
#SBATCH --output=outFiles/2.1metatreeB.%A.out # File to which standard out will be written. %A is replaced by the job ID
#SBATCH --array=1-1                           # Job array range

# Print time stamp and other info for output files
date;hostname;pwd

echo "Set variables"
echo "Make the arguments for the sbatch command in the lower loop"
# These arguments come from 1_pplacer_A.sh script
GENE=$1
NCORES=$2
CUTOFF=$3
DATASET=$4
GENE_DIR=$5
# Variables that take to the where the metaT sample files and what samples to work on
SAMP_LIST="data/sampleLists/${DATASET}.sampleList.txt"
SAMP_DIR="data/${DATASET}"

echo "Loop over all samples listed in SAMP_LIST"
# We are calling 1_pplacer_C.${DATASET}.sh and passing the 1_pplacer_A.sh arguments to it  
for subject in $(cat ${SAMP_LIST}); do
  sbatch 1_pplacer_C_${DATASET}.sh $GENE $NCORES $CUTOFF $DATASET $GENE_DIR $SAMP_LIST $SAMP_DIR $subject
done