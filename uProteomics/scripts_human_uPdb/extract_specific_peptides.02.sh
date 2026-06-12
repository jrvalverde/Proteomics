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
# HEAVILY DOCUMENTED AND POSSIBLY, EVEN KEPT IN-HOUSE.
#
#   See version 03 of both.
#

db=${1:-db_pep_8-30.1}

for i in $db/*.faa ; do
    fa=${i##*/}
    cl=${db}/${fa}.clstr
    mkdir -p ${db}.sp
    ou=${db}.sp/${fa}  # the final output fasta file
    echo $fa $cl
    if [ -s /${db}.sp/${fa} ] ; then continue ; fi  # already done
    # we need a specific output directory so we do not
    # overwrite each other
    truncate -s 0 $ou       # create an empty ouput fasta file
    
    perl ./extract_Nseq_clstrs.pl $i $cl tmp/$fa
    # We exploit here the default cluster size setting (4th argument), 
    # which is 1.
    # This will generate one fasta file per cluster in tmp
    # named with the cluster number (0 1 2 3 ... 100 101 ...)
    # This is so because it first builds a list of clusters,
    # and then runs over the FASTA file and saves each sequence
    # to its own cluster FASTA file, so it needs to make a separate
    # FASTA file for each cluster as cluster FASTAs will be appended
    # incrementally as their sequences are found in the original FASTA.
    #
    #   (We could modify it to save to only one file is size_cutoff is 1)
    #
    # We need to concatenate them now
    #/bin/ls -v tmp/$fa | while read name ; do 
    find tmp/$fa -print | sort -n | while read name ; do 
        #cho $name
        cat $name >> $ou    # append cluster fasta to output fasta file
        rm $name
    done
    rm -rf ./tmp/${fa}
    echo $i processed
done
