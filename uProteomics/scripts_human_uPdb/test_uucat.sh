cd uuprots

if [ ! -s ../uucat/all.fa ] ; then
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
        > ../uucat/all.faa 
        # all.faa   : 12 954 476
        # all.1.faa :  9 274 434
fi

if [ ! -s ../uucat/openprot.faa ] ; then
    cat openprot* \
        > ../uucat/openprot.faa
    # openprot.faa   : 134 009
    # openprot.1.faa :  55 256

fi

if [ ! -s ../uucat/openp+unip.faa ] ; then
    cat openprot* \
        uniprot* \
        > ../uucat/openp+unip.faa
    # openp+unip.faa   : 138 128
    # openp+unip.1.faa :  55 469
fi

if [ ! -s ../uucat/ncbi+openp+unip.faa ] ; then
    cat ncbi* \
        openprot* \
        uniprot* \
        > ../uucat/ncbi+openp+unip.faa
    # cbi+openp+unip.faa   : 138 151
    # ncbi+openp+unip.1.faa :  55 470
fi

if [ ! -s ../uucat/ncbi+ncEP+openp+unip.faa ] ; then
    cat nc* \
        openprot* \
        uniprot* \
        > ../uucat/ncbi+ncEP+openp+unip.faa
    # ncbi+ncEP+openp+unip.faa   : 138 154
    # ncbi+ncEP+openp+unip.1.faa :  55 471
fi

if [ ! -s ../uucat/ncbi+ncEP+openp+smp+unip.faa ] ; then
    cat nc* \
        openprot* \
        SmProt* \
        uniprot* \
        > ../uucat/ncbi+ncEP+openp+smp+unip.faa 
    # ncbi+ncEP+openp+smp+unip.faa   : 478 879
    # ncbi+ncEP+openp+smp+unip.1.faa : 139 296
fi

if [ ! -s ../uucat/ncbi+ncEP+openp+smp+spencer+unip.faa ] ; then
    cat nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../uucat/ncbi+ncEP+openp+smp+spencer+unip.faa
    # ncbi+ncEP+openp+smp+spencer+unip.faa   : 531 887
    # ncbi+ncEP+openp+smp+spencer+unip.1.faa : 143 439
fi

if [ ! -s ../uucat/ensembl+ncbi+ncEP+openp+unip.faa ] ; then
    cat ensembl* \
        nc* \
        openprot* \
        uniprot* \
        > ../uucat/ensembl+ncbi+ncEP+openp+unip.faa
    # ensembl+ncbi+ncEP+openp+unip.faa    :140 230
    # ensembl+ncbi+ncEP+openp+unip.1.faa :  55 628
fi

if [ ! -s ../uucat/ensembl+ncbi+openp+smp+spencer+unip.faa ] ; then
    cat ensembl* \
        nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../uucat/ensembl+ncbi+openp+smp+spencer+unip.faa
    # ensembl+ncbi+openp+smp+spencer+unip.faa   : 533 963
    # ensembl+ncbi+openp+smp+spencer+unip.1.faa : 143 592

fi


if [ ! -s ../uucat/cnc+ensembl+ncbi+openp+smp+spencer+unip.faa ] ; then
    cat cnc* \
        ensembl* \
        nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../uucat/cnc+ensembl+ncbi+openp+smp+spencer+unip.faa
    # cnc+ensembl+ncbi+openp+smp+spencer+unip.faa   : 535 715
    # cnc+ensembl+ncbi+openp+smp+spencer+unip.1.faa : 144 152
fi

if [ ! -s ../uucat/cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa ] ; then
    cat cnc* \
        ensembl* \
        MicroProteinDB* \
        nc* \
        openprot* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../uucat/cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa
    # cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa   : 571 349
    # cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.1.faa : 167 358
fi


# So far, these are relatively small, let us check now without openprot

if [ ! -s ../uucat/cnc+ensembl+micro+ncbi+smp+spencer+unip.faa ] ; then
    cat cnc* \
        ensembl* \
        MicroProteinDB* \
        nc* \
        SmProt* \
        SPENCER* \
        uniprot* \
        > ../uucat/cnc+ensembl+micro+ncbi+smp+spencer+unip.faa
    # cnc+ensembl+micro+ncbi+smp+spencer+unip.faa   : 437 340
    # cnc+ensembl+micro+ncbi+smp+spencer+unip.1.faa : 119 013
fi


cd ../uucat

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

