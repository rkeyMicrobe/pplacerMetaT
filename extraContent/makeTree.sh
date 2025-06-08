#!/bin/bash

#SBATCH --job-name=outFiles/makeTrees.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=8GB
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH -t 48:00:00
#SBATCH --partition=hpg-milan
#SBATCH --output=outFiles/%a.makeTrees.%A.out
#SBATCH --array=1-2

# Print time stamp and other info for output files
date
cd /blue/b.durham/rebeccakey/4.1_treeMaker
hostname;pwd

echo "Set variables"
# Variables
GENE=$(sed -n ${SLURM_ARRAY_TASK_ID}p geneList.txt)
GENE_PATH="genes/${GENE}"

# ----------------------------------------------------------
# ----------------------------------------------------------

# PREPARE YOUR GENE CANDIDATES (you need two starting files: ingroup and outgroup txt files per gene)
echo "WORKING ON ${GENE}"

# Find orthologs IDs to your ingroup queries in curated databases
module purge; ml ncbi_blast/2.10.1
echo "Blast reference databases with gene ingroup file"
date

blastp -query "${GENE_PATH}/start/${GENE}_query_ingroup.txt" \
-db tools/db/MARINEREFII/blastMarRefII/MarRefIIaa \
-out "${GENE_PATH}/${GENE}_blast_marDB.txt" \
-evalue 1e-10 \
-outfmt "6 qseqid sseqid evalue stitle" \
-max_target_seqs 1000

blastp -query "${GENE_PATH}/start/${GENE}_query_ingroup.txt" \
-db tools/db/opistho_blastdb/opistho_aa \
-out "${GENE_PATH}/${GENE}_blast_opisDB.txt" \
-evalue 1e-10 \
-outfmt "6 qseqid sseqid evalue stitle"

# Pull the sequences and classification info that correspond to the IDs
module purge; ml python/2.7.18
echo "Get fasta sequences from each database and name them"
date

echo "doing the MarRef..."
python tools/create_fasta.py "${GENE_PATH}/${GENE}_blast_marDB.txt" "${GENE_PATH}/${GENE}_fasta_marDB.txt"
python tools/db/MARINEREFII/Rename_MarineRefII_Seqs.py \
-I tools/db/MARINEREFII/MarineRefII_seqIDinfo.csv \
-S "${GENE_PATH}/${GENE}_fasta_marDB.txt"
mv "named_${GENE}_fasta_marDB.txt" "${GENE_PATH}"

echo "doing the Opisthokont..."
python tools/create_fasta_opistho.py "${GENE_PATH}/${GENE}_blast_opisDB.txt" "${GENE_PATH}/${GENE}_fasta_opisDB.txt"
python tools/fnames2.py \
-t "${GENE_PATH}/${GENE}_fasta_opisDB.txt"

# Merge blast hits to one file; Merge your starting query seqs (both in and outgroup) to one file
echo "Concatenating fasta sequences"
date
cat "${GENE_PATH}/named_${GENE}_fasta_marDB.txt" "${GENE_PATH}/${GENE}_fasta_opisDB.fn.txt" > "${GENE_PATH}/${GENE}_fasta.txt"

echo "Appending outgroup sequences onto final"
cat "${GENE_PATH}/start/${GENE}_query_outgroup.txt" "${GENE_PATH}/start/${GENE}_query_ingroup.txt" > "${GENE_PATH}/${GENE}_query_refs.txt"

# Clean up your sequences so all are the same length and structure
echo "Removing sequences by length"
python tools/sequence_cleaner_rsk.py "${GENE_PATH}/${GENE}_fasta.txt" 100 1500

echo "Getting rid of unnecessary chars"
grep -o '>' "${GENE_PATH}/${GENE}_fasta.txt" | wc -l
grep -o '>' "${GENE_PATH}/${GENE}_fasta_clean.txt" | wc -l

# Merge blast hits and your queries to one master file for analysis
echo "Merging both your queries and blasted sequences into one fasta file"
cat "${GENE_PATH}/${GENE}_query_refs.txt" "${GENE_PATH}/${GENE}_fasta_clean.txt" > "${GENE_PATH}/${GENE}_all.fasta"
date


# ----------------------------------------------------------
# ----------------------------------------------------------

# PREPARE SEQUENCE ALIGNMENTS FOR RAXML ANALYSIS

# Cluster your sequences to reduce computation time and alignment errors
module purge; ml usearch/11.0.667
echo "Using usearch to cluster similar sequences"
date
usearch -cluster_smallmem "${GENE_PATH}/${GENE}_all.fasta" \
-id 0.8 -sortedby other -centroids "${GENE_PATH}/${GENE}_clust80.fasta"

usearch -cluster_smallmem "${GENE_PATH}/${GENE}_all.fasta" \
-id 0.8 -sortedby other -uc "${GENE_PATH}/${GENE}_clust80.uc"

# Perform your alignments
module purge; ml mafft/7.505
echo "Using mafft to perform alignments"
date
mafft --thread 8 --genafpair --maxiterate 16 --leavegappyregion \
--inputorder "${GENE_PATH}/${GENE}_clust80.fasta" > "${GENE_PATH}/${GENE}_mafft.fasta"

# Trim your alignments
module purge; ml trimal/1.4.1
echo "Using trimAl, trim the alignments"
date
trimal -in "${GENE_PATH}/${GENE}_mafft.fasta" -out "${GENE_PATH}/${GENE}_trimal1.fasta" -gt .05
trimal -in "${GENE_PATH}/${GENE}_trimal1.fasta" -out "${GENE_PATH}/${GENE}.trimal2.fasta" -resoverlap 0.5 -seqoverlap 50

echo "Remove trailing length tag on trimal2.fasta file"
head "${GENE_PATH}/${GENE}.trimal2.fasta"
sed -i 's/ 822 bp//g' "${GENE_PATH}/${GENE}.trimal2.fasta"

# Convert your fasta to a useable input for prottest3 and raxml commands
echo "Converting fasta to phylip format"
date
perl tools/Fasta2Phylip.pl "${GENE_PATH}/${GENE}.trimal2.fasta" "${GENE_PATH}/${GENE}_phylip.phy"

# Get a best-fit substitution model for your alignments
# Note: The substitution model describes the rates at which amino acid substitutions occur during evolution and takes into account factors such as amino acid properties, structural constraints, and evolutionary processes
module purge; ml prottest3/3.4.2
echo "Finding the best model with prottest"
date
prottest3 -i "${GENE_PATH}/${GENE}_phylip.phy" -S 0 -JTT -LG -DCMut -Dayhoff -WAG -Blosum62 -VT -IG -F -BIC -threads 8 -o "${GENE_PATH}/${GENE}_prottest.out"

# ----------------------------------------------------------
# ----------------------------------------------------------

# MAKE THE TREES USING RAXML

# Feed the result alignments into the raxml command 
module purge; ml raxml
echo "Make a fast tree with raxml"
date
cd ${GENE_PATH}
pwd
raxmlHPC-PTHREADS-SSE3 -T 8 -f E -p 287 -m "PROTGAMMAILG" -n "${GENE}.tre" -s "${GENE}_phylip.phy" -O
raxmlHPC-PTHREADS-SSE3 -T 8 -f e -m "PROTGAMMAILG" -t RAxML_fastTree."${GENE}.tre" -n "brL.${GENE}.tre" -s "${GENE}_phylip.phy" -O

# Push the input files needed for pplacer analysis using ocean water transcripts
echo "Push result outputs to a fresh directory for pplacer analysis"
date
cp "${GENE}.trimal2.fasta" "../../trees/${GENE}"
cp "RAxML_info.brL.${GENE}.tre" "../../trees/${GENE}"
cp "RAxML_result.brL.${GENE}.tre" "../../trees/${GENE}"

cd /blue/b.durham/rebeccakey/4.1_treeMaker
echo "END"
date; pwd 


# dead code
# Reciprocal test against Thps
#blastp -query "${GENE}.cat.fasta.txt" -db ~/Desktop/bpdurham/Thaps3_blast/blastdb/Thaps_aa -out "${GENE}.recip.txt" -evalue 1 -outfmt "6 qseqid sseqid evalue stitle" -max_target_seqs 1

