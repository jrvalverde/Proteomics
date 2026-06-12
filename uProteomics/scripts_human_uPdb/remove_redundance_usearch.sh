#!/bin/bash
#
# This is much-much-much faster than clustering with cd-hit
#   and apparently gives the same (or almost the same) results 
#   (although we also request to include cluster-size in the 
#   annotation)
#

db=${1}

mkdir -p ${db}.1
cd $db

for i in ./*.faa ; do
    echo "Deduplicating $i"
    # do not repeat work already done or in progress
    if [ -s ../${db}.1/$i ] ; then continue ; fi
    date
    time ../bin/usearch12 -fastx_uniques $i \
        -fastaout ../${db}.1/${i} \
        -sizeout
    date
done |& tee ../${db}.1/uniq-usearch12.log
cd -

exit

# the fasta output files have all the centroids and their corresponding
# cluster size, if we want a fasta file containing only unique sequences
# we need to filter sequences annotated with ;size=1;
#   note that we expect all-cluster files to be named .fasta
cd ${db}.1
for i in *.fasta ; do 
    echo $i ;
    cat $i \
    | tr '\n' '@' \
    | sed -e 's/@>/\n>/g' \
    | grep ';size=1;' \
    | tr '@' '\n' \
    > ${i%sta}a ; 
done

