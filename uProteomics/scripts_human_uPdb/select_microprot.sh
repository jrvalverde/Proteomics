#!/bin/bash
#
# to be run in the data directory
#

topdir=`pwd`
fa_db=$topdir/fa.db
fa_uprot_db=$topdir/fa.uprot.db
uprot_db=$fa_uprot_db/uprots

mkdir -p uprots
#rm -f uprots/*

min=8
max=100
# should be max=101

mkdir -p uprots

for db in cnc* ens* Lnc* Met* Mi* nc* open* sm* Sm* SP* uni* uorf* ; do
    if [ ! -d $db/*tagged ] ; then continue ; fi

    cd $db/*tagged
    pwd

    for i in *.fa* ; do
        name=${i%.*}
        if [ ! -f $i ] ; then continue ; fi
        out=$uprot_db/${db}-${name}.${min}-${max}.faa
        if [ ! -s "$out" ] ; then
  	        seqkit seq -g -m $min -M $max "$i" > "$out"
            echo "$i"
	        echo -n "    Original: "
	        grep -c '>' $i
	        echo -n "    Filtered: "
	        grep -c '>' $out
            echo ""
        fi
    done

    cd -

done

# This has selected all microproteins from the databases
#
# Next repeat in all proteins in the db directory

for db in $fa_db ; do

    cd $db
    pwd

    mkdir -p ../uprots

    for i in *.fa* ; do
        name=${i%.*}
        out=$uprot_db/${db}-${name}.${min}-${max}.faa
        if [ ! -f $i ] ; then continue ; fi
        if [ -s $out ] ; then continue ; fi
        if [ ! -s "$out" ] ; then
            echo "Processing $i to $out"
  	        seqkit seq -g -m $min -M $max "$i" > "$out"
            echo "$i"
	        echo -n "    Original: "
	        grep -c '>' $i
	        echo -n "    Filtered: "
	        grep -c '>' $out
            echo ""
        else
            echo "$i already done"
        fi
    done

    cd -

done

