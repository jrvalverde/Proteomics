cd uprots

if [ ! -s ../ucat/all.fa ] ; then
    cat cncRNADB* \
        ensembl* \
        LncPep* \
        MetamORF* \
        MicroProteinDB* \
        MiTPeptideDB* \
        ncbi* \
        ncEP* \
        openprot* \
        smORFunction* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../ucat/all.faa 
        # all.faa   : 10 798 045
        # all.1.faa : 19 107 111
fi

if [ ! -s ../ucat/openprot.faa ] ; then
    cat openprot* \
        > ../ucat/openprot.faa
    # openprot.faa   : 1 620 162
    # openprot.1.faa :   599 880

fi

if [ ! -s ../ucat/openp+unip.faa ] ; then
    cat openprot* \
        uniprot* \
        > ../ucat/openp+unip.faa
    # openp+unip.faa   : 1 660 395
    # openp+unip.1.faa :   601 020
fi

if [ ! -s ../ucat/ncbi+openp+unip.faa ] ; then
    cat ncbi* \
        openprot* \
        uniprot* \
        > ../ucat/ncbi+openp+unip.faa
    # cbi+openp+unip.faa   : 1 663 202
    # ncbi+openp+unip.1.faa :   601 062
fi

if [ ! -s ../ucat/ncbi+ncEP+openp+unip.faa ] ; then
    cat nc* \
        openprot* \
        uniprot* \
        > ../ucat/ncbi+ncEP+openp+unip.faa
    # ncbi+ncEP+openp+unip.faa   : 1 663 216
    # ncbi+ncEP+openp+unip.1.faa :   601 064
fi

if [ ! -s ../ucat/ncbi+ncEP+openp+smp+unip.faa ] ; then
    cat nc* \
        openprot* \
        SmProt* \
        uniprot* \
        > ../ucat/ncbi+ncEP+openp+smp+unip.faa 
    # ncbi+ncEP+openp+smp+unip.faa   : 2 554 371
    # ncbi+ncEP+openp+smp+unip.1.faa :   716 773
fi

if [ ! -s ../ucat/ncbi+ncEP+openp+smp+spencer+unip.faa ] ; then
    cat nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../ucat/ncbi+ncEP+openp+smp+spencer+unip.faa
    # ncbi+ncEP+openp+smp+spencer+unip.faa   : 2 613 348
    # ncbi+ncEP+openp+smp+spencer+unip.1.faa :   719 508
fi

if [ ! -s ../ucat/ensembl+ncbi+ncEP+openp+unip.faa ] ; then
    cat ensembl* \
        nc* \
        openprot* \
        uniprot* \
        > ../ucat/ensembl+ncbi+ncEP+openp+unip.faa
    # ensembl+ncbi+ncEP+openp+unip.faa   : 1 696 332
    # ensembl+ncbi+ncEP+openp+unip.1.faa :   608 359
fi

if [ ! -s ../ucat/ensembl+ncbi+openp+smp+spencer+unip.faa ] ; then
    cat ensembl* \
        nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../ucat/ensembl+ncbi+openp+smp+spencer+unip.faa
    # ensembl+ncbi+openp+smp+spencer+unip.faa   : 2 646 464
    # ensembl+ncbi+openp+smp+spencer+unip.1.faa :   726 736

fi


if [ ! -s ../ucat/cnc+ensembl+ncbi+openp+smp+spencer+unip.faa ] ; then
    cat cnc* \
        ensembl* \
        nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../ucat/cnc+ensembl+ncbi+openp+smp+spencer+unip.faa
    # cnc+ensembl+ncbi+openp+smp+spencer+unip.faa   : 2 648 235
    # cnc+ensembl+ncbi+openp+smp+spencer+unip.1.faa :   727 143
fi

if [ ! -s ../ucat/cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa ] ; then
    cat cnc* \
        ensembl* \
        MicroProteinDB* \
        nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../ucat/cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa
    # cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa   : 2 698 188
    # cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.1.faa :   746 983
fi

# So far, these are relatively small, let us check now without openprot

if [ ! -s ../ucat/cnc+ensembl+micro+ncbi+smp+spencer+unip.faa ] ; then
    cat cnc* \
        ensembl* \
        MicroProteinDB* \
        nc* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../ucat/cnc+ensembl+micro+ncbi+smp+spencer+unip.faa
    # cnc+ensembl+micro+ncbi+smp+spencer+unip.faa   : 1 078 026
    # cnc+ensembl+micro+ncbi+smp+spencer+unip.1.faa :	238 651
fi




cd ../ucat

for i in *[A-Za-z].faa ; do
    name=${i%.faa}
    echo $name
    #continue
    if [ ! -s ${name}.1.faa ] ; then
	file=$i
	name=${file%.faa}
	ext=${file##*.}

	cd-hit -i $file -o $name.1.$ext -c 1 -M 100000 -T 24 -n 4 -l 8
    fi
done

