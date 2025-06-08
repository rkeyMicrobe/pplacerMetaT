GENE=$1

CUTOFF="40"
NCORES="12"
SUBJECT_DIR="/mnt/nfs/ryan/Gradients1/mORFeus_v2"
SUBJECT_LIST="/mnt/nfs/ryan/Gradients1/mORFeus_v2/gradients1_morfeus_handles.txt"
TAX_DB="/mnt/nfs/home/bpdurham/taxonomy.db"

#/mnt/nfs/home/bpdurham/tree_commands/make_seq_info_all.py $GENE.trimal2.fasta
#/mnt/nfs/home/bpdurham/taxtastic-0.5.4/taxit update_taxids -d $TAX_DB -o seq_info_all.updated.csv seq_info_all.csv
#cat seq_info_all.updated.csv | awk -F, '{print $2}' | sort | uniq | sed '/tax_id/d' > tax_ids_all.txt
#/mnt/nfs/home/bpdurham/taxtastic-0.5.4/taxit taxtable -d $TAX_DB -t tax_ids_all.txt -o taxa.csv

#make the refpkg with this:
#/mnt/nfs/home/bpdurham/taxtastic-0.5.4/taxit create -l $GENE -P $GENE.refpkg \
#--taxonomy taxa.csv \
#--aln-fasta $GENE.trimal2.fasta \
#--seq-info seq_info_all.updated.csv \
#--tree-stats RAxML_info.brL.$GENE.tre \
#--tree-file RAxML_result.brL.$GENE.tre \
#--no-reroot

#seqmagick convert $GENE.trimal2.fasta $GENE.sto
#hmmbuild $GENE.hmm $GENE.sto

#mkdir pplacer_gradients1; cd pplacer_gradients1/

REFPKG="../$GENE.refpkg"

/mnt/nfs/home/bpdurham/tree_commands/pplacer_script_gradients.sh $GENE $NCORES $CUTOFF

tar czf $GENE.gradients.jplace.tar.gz *jplace
mkdir jplace; mv *jplace jplace; cd jplace/

# make the 'all' tree:
guppy fat -o $GENE.all.xml --average $GENE.*.jplace
/mnt/nfs/home/bpdurham/tree_commands/treecolor2.py -o $GENE.allfatcolor.xml $GENE.all.xml 
/mnt/nfs/home/bpdurham/tree_commands/treecolor2.py -o $GENE.allcolor.xml -a $GENE.all.xml 
mv *all* ..

# getting a sample list:
ls $GENE.*jplace | awk -F. {'print $2'} | sort | uniq > ../sample_list.txt

# process CSVs by summed samples:
for sample in $(cat ../sample_list.txt); do
  ~/Desktop/bpdurham/pplacer/guppy to_csv -o $GENE.$sample.taxID.csv $GENE.$sample.*.jplace
done
mkdir ../summed_csv; mv *csv ../summed_csv; cd ../summed_csv




#TREECOLORS="/Users/bryndandurham/Desktop/bpdurham/tree_commands/treecolor_files/treecolors_w_proks.csv"
#ingroup_list="$GENE.just_edges.txt"
#NORMFACTORS="/Users/bryndandurham/Desktop/bpdurham/tree_commands/gradients1.norm_factor_SUMS.csv"

#~/Desktop/bpdurham/tree_commands/count_pplacer_csv_by_taxonomy.py -e -g -c $TREECOLORS -n $NORMFACTORS $(ls $GENE.*.taxID.csv | head -1) > $GENE.counts_results_gradients1.csv
#for csv in $(ls $GENE.*.taxID.csv); do ~/Desktop/bpdurham/tree_commands/count_pplacer_csv_by_taxonomy.py -i $ingroup_list -g -c $TREECOLORS -n $NORMFACTORS $csv >> $GENE.counts_results_gradients1.csv; done




#get the tree-edge structure from any of the jplace files:
#head -2 $GENE.*.jplace > $GENE.tree.txt

#manually excise the desired structure and paste into a file (e.g. $GENE.raw_edges.txt)

#Pull out edge numbers from the 'cut' tree (this seems to work on Mac but not bloom!!)
#grep -o -E "\{\d+\}" $GENE.raw_edges.txt | tr -d '{}' > $GENE.just_edges.txt