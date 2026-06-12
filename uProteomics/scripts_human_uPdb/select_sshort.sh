#!/bin/bash
#
# to be run in the data directory
#

mkdir -p uuprots
#rm -f uuprots/*

here=`pwd`
min=8
max=30

for db in cnc* ens* Lnc* Met* Mi* nc* open* sm* Sm* SP* uni* uorf* ; do
    if [ ! -d $db/*sados ] ; then continue ; fi

    cd $db/*sados
    pwd

    mkdir -p ../uuprots

    for i in *.fa* ; do
        name=${i%.*}
        if [ ! -f $i ] ; then continue ; fi
        out=${db}-${name}.${min}-${max}.faa
	    seqkit seq -g -m $min -M $max $i > $here/uuprots/$out
        echo "$i"
	    echo -n "    Original: "
	    grep -c '>' $i
	    echo -n "    Filtered: "
	    grep -c '>' $here/uuprots/$out
        echo ""
    done

    cd -

done
    
    
