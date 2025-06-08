HANDLE=$1
NCORES=$2 # e.g., 16 for gross and 8 for match
BITSCORE_CUTOFF=$3 # original default value was 35

SUBJECT_LIST="/mnt/nfs/ryan/Gradients1/mORFeus_v2/gradients1_morfeus_handles.txt"
SUBJECT_DIR="/mnt/nfs/ryan/Gradients1/mORFeus_v2"

for subject in $(cat $SUBJECT_LIST); do
  # build the filepath for subject files
  SUBJECT_FILE="$subject".6tr.orfs40.fasta.gz
  # do the hmmsearch on each
  hmmsearch -A "$HANDLE"."$subject".query.sto -T "$BITSCORE_CUTOFF" --incT "$BITSCORE_CUTOFF" --cpu "$NCORES" --tblout "$HANDLE".hmm_out.tab ../"$HANDLE".hmm $SUBJECT_DIR/$SUBJECT_FILE
  # output: "$HANDLE"."$subject".query.sto
  # options for modification: change E-value threshold
  # -E <x>  where x > 0
  # use parallel processing: --cpu 16

  # use hmmalign to align query hits to the reference alignment
  hmmalign -o "$HANDLE"."$subject".aln.sto --mapali ../"$HANDLE".sto ../"$HANDLE".hmm "$HANDLE"."$subject".query.sto
  # output: "$HANDLE"."$subject".aln.sto
  
  seqmagick convert "$HANDLE"."$subject".aln.sto "$HANDLE"."$subject".aln.fasta
  
  # run pplacer using refpkg
  pplacer -c ../"$HANDLE".refpkg --keep-at-most 1 "$HANDLE"."$subject".aln.fasta
  # output:

done
