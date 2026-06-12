#!/bin/bash
#
#   For generic clusters we need the original fasta, because the
# clustered fasta only contains the representative sequence, but
# not the other cluster sequences, but for specific peptides 
# (peptides that belong only to one protein and therefore form 
# size-1 clusters) we can save time using the clustered fasta because
# we only need the centroid/representative in it.
#
#   We can also save time if we prepare things so that we can run 
# many scripts in parallel in separate machines.
#
#   We can further improve speed by changing extract_Nseq_clstrs.pl
# so that when size_cutoff is 1, the output destination is taken 
# not as a directory to save clusters separately, but as a file to
# save all size-1 clusters together. THIS IS INCONSISTENT BEHAVIOUR
# COUNTERINTUITIVE AND DANGEROUS, SO THE MODIFIED SCRIPT SHOULD BE
# HEAVILY DOCUMENTED, BUT THE SPEEDUP GAINED CAN BE IMPRESSIVE.
#
#   IN THIS VERSION WE DO USE THE MODIFIED extract_Nseq_clstrs.03.pl
#

db=${1:-db_pep_8-30.1}

for i in $db/*.fasta ; do
    fa=${i##*/}
    cl=${i}.clstr
    ou=${i%sta}a  # the final output fasta file
    echo $fa $cl
    if [ -s ${ou} ] ; then continue ; fi  # already done
    
    # ensure the output file exists to avoid race conditions
    # when run in parallel
    truncate -s 0 $ou       # create an empty ouput fasta file
    
    perl ./extract_Nseq_clstrs.03.pl $i $cl ${ou} 1
    # We exploit here the cluster size setting (4th argument), 
    # which is set to 1.
    # This makes use of a specific version of the perl script that
    # distinguishes when size_cutoff is one, and then, instead of
    # using separate files for each cluster, saves all size-1 clusters
    # to a single file
done
