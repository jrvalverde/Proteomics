#!/bin/bash

in=${1:-ucat.1_pep_8-30.1}       # directory with peptide databases to analyze

sc=${in}_scored
de=${in}_detectable

mkdir -p $sc
mkdir -p $de

for i in $in/mabad*.faa ; do
    name=${i##*/}
  # we'll do it in parallel because it is mostly computationally intensive
  (
    if [ ! -s $i ] ; then continue ; fi
    echo "Processing $name"
    # predict peptide detectability and add to annotation comments
    if [ ! -s $sc/$name ] ; then
        bin/detectability/predict_detectability.sh $i $sc/$name
    fi
    if [ ! -s $de/$name ] ; then
        # convert to one single line substituting '\n' by '@'
        # insert a '\n' before '>' effectively having each sequence in one line
        # filter lines (sequences) corresponding to observable peptides
        # put back the '\n' instead of '@'
        # and save in detectable directory
        cat $sc/$name \
        | tr '\n' '@' \
        | sed -e 's/@>/\n>/g' \
        | grep 'class=MObs' \
        | tr '@' '\n' \
        > $de/$name
    fi
  ) >& log/${name}_detectability.log &
done
wait

