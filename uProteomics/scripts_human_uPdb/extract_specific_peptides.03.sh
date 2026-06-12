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

for i in $db/*.faa ; do
    fa=${i##*/}
    cl=${db}/${fa}.clstr
    mkdir -p ${db}.sp
    ou=${db}.sp/${fa}  # the final output fasta file
    echo $fa $cl
    if [ -s ${db}.sp/${fa} ] ; then continue ; fi  # already done
    
    # ensure the output file exists to avoid race conditions
    truncate -s 0 $ou       # create an empty ouput fasta file
    
    if [ ! -s $cl ] ; then
        # we must have used usearch12 which currently cannot produce
        # an uc file, but then the fasta must have an annotated
        # cluster size:
        # the fasta output files have all the centroids and their 
        # corresponding cluster size, if we want a fasta file 
        # containing only unique sequences we need to filter sequences 
        # annotated with ;size=1;
        #   convert '\n' to '@' and make a one line file
        #   change '@>' to '\n>' and put each entry@sequence in one line
        #   select lines corresponding to size=1 clusters
        #   restore all '\n' and save
        echo $i ;
        cat $i \
        | tr '\n' '@' \
        | sed -e 's/@>/\n>/g' \
        | grep ';size=1;' \
        | tr '@' '\n' \
        > $ou ; 
    else
        perl ./extract_Nseq_clstrs.03.pl $i $cl ${ou}
        # We exploit here the default cluster size setting (4th argument), 
        # which is 1.
        # This makes use of a specific version of the perl script that
        # distinguishes when size_cutoff is one, and then, instead of
        # using separate files for each cluster, saves all size-1 clusters
        # to a single file
    fi
done
