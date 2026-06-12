#!/bin/bash
#
#   TO BE RUN IN THE TOP DIRECTORY ABOVE THE fa.*.db DIRECTORIES


topdir=`pwd`
fa_uprot_db=$topdir/fa.uprot.db
fa_uprot_multi_db=$topdir/fa.uprot.multi.db

uprot_db_1=$fa_uprot_db/uprots.1      # unique microproteins (non redundant)
ucat_db=$fa_uprot_multi_db/ucat
ucat_db_1=$fa_uprot_multi_db/ucat.1

mkdir -p $ucat_db
mkdir -p $ucat_db_1

cd $uprot_db_1


# These are only proteins < 100 a.a. (after removing redundance)
#   There may be a problem with proteins containing a '*' in their
#   sequence, but it should be automatically fixed in later steps
#

mfa=all.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat cncRNADB*.faa \
        ensembl*.faa \
        LncPep*.faa \
        MetamORF*.faa \
        MicroProteinDB*.faa \
        MiTPeptideDB*.faa \
        ncbi*.faa \
        ncEP*.faa \
        openprot*.faa \
        smORFunction*.faa \
        SmProt*.faa \
        SPENCER*.faa \
        uniprot*.faa \
        db*.faa \
        > $ucat_db/$mfa 
fi

mfa=openprot.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat openprot*.faa \
        > $ucat_db/$mfa
fi

mfa=openp+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat openprot*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=ncbi+openp+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat ncbi*.faa \
        openprot*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=ncbi+ncEP+openp+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat nc*.faa \
        openprot*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=ncbi+ncEP+openp+smp+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat nc*.faa \
        openprot*.faa \
        SmProt*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=ncbi+ncEP+openp+smp+spencer+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat nc*.faa \
        openprot*.faa \
        SmProt*.faa \
        SPENCER*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=ensembl+ncbi+ncEP+openp+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat ensembl*.faa \
        nc*.faa \
        openprot*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=ensembl+ncbi+openp+smp+spencer+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat ensembl*.faa \
        nc*.faa \
        openprot*.faa \
        SmProt*.faa \
        SPENCER*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=cnc+ensembl+ncbi+openp+smp+spencer+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat cnc*.faa \
        ensembl*.faa \
        nc*.faa \
        openprot*.faa \
        SmProt*.faa \
        SPENCER*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=cnc+ensembl+micro+ncbi+openp+smp+spencer+unip.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat cnc*.faa \
        ensembl*.faa \
        MicroProteinDB*.faa \
        nc*.faa \
        openprot*.faa \
        SmProt*.faa \
        SPENCER*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=RNAseq.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-RNASeq_Uniprot_human*.faa \
        db-human_Uniprot_v3_RNASeq*.faa \
        cncRNADB-Translated_ncRNA*.faa \
        SmProt-SmProt2_circRNAtranslation*.faa \
        > $ucat_db/$mfa
fi

mfa=small.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        ensembl-Homo_sapiens.GRCh38.pep.all*.faa \
        MicroProteinDB-sequence*.faa \
        ncbi-refseq-human*.faa \
        ncEP-*.faa \
        SmProt-SmProt_*.faa \
        SPENCER-*.faa \
        uniprot-*.faa \
        > $ucat_db/$mfa
fi

mfa=small2.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        ensembl-Homo_sapiens.GRCh38.pep.all*.faa \
        MicroProteinDB-sequence*.faa \
        ncbi-refseq-human*.faa \
        ncEP-*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_MS*.faa \
        SPENCER-*.faa \
        uniprot-*.faa \
        > $ucat_db/$mfa
fi

mfa=small2+naug.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        ensembl-Homo_sapiens.GRCh38.pep.all*.faa \
        MicroProteinDB-sequence*.faa \
        ncbi-refseq-human*.faa \
        ncEP-*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        uniprot-*.faa \
        > $ucat_db/$mfa
fi

mfa=small2+ribo.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        ensembl-Homo_sapiens.GRCh38.pep.all*.faa \
        MicroProteinDB-sequence*.faa \
        ncbi-refseq-human*.faa \
        ncEP-*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo* \
        SPENCER-*.faa \
        uniprot-*.faa \
        > $ucat_db/$mfa
fi

mfa=small2+naug+ribo.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        ensembl-Homo_sapiens.GRCh38.pep.all*.faa \
        MicroProteinDB-sequence*.faa \
        ncbi-refseq-human*.faa \
        ncEP-*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Ribo*.faa\
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        uniprot-*.faa \
        > $ucat_db/$mfa
fi

mfa=mprot+smp2+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mprot+smp2+naug+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mprot+smp2+ribo+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo* \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mprot+smp2+naug+ribo+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+naug+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+ribo+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo* \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=nc+uniprot.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat nc*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=cnc+nc+uniprot.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat cnc*.faa \
        nc*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=cnio.faa
# we avoid db-RNASeq_Uniprot_human.8-100.faa because it is too large (520831)
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-human_Uniprot_v3_RNASeq.8-100.faa \
        db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        > $ucat_db/$mfa
fi

mfa=cnio+.faa
# we avoid db-RNASeq_Uniprot_human.8-100.faa because it is too large (520831)
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-*.faa \
        > $ucat_db/$mfa
fi

mfa=cnio+nc.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-human_Uniprot_v3_RNASeq.8-100.faa \
        db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        nc*.faa \
        > $ucat_db/$mfa
fi

mfa=cnio+uniprot.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-human_Uniprot_v3_RNASeq.8-100.faa \
        db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi

mfa=cnio+nc+uniprot.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-human_Uniprot_v3_RNASeq.8-100.faa \
        db-hUniprot_20373_sorf_*.faa \
        db-hUniprot_23candidates*.faa \
        db-smprot_openprot_SEPs_MAbad_c*.faa \
        nc*.faa \
        uniprot*.faa \
        > $ucat_db/$mfa
fi
   
mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+cnc+oprot2pepp13.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        cncRNADB-Translated_ncRNA.8-100.faa \
        openprot-human-openprot-2_2-refprots+altprots+isoforms-min_2_pep-GRCh38.p13-uniprot2022_06_01.8-100.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+oprot2pepp13.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        openprot-human-openprot-2_2-refprots+altprots+isoforms-min_2_pep-GRCh38.p13-uniprot2022_06_01.8-100.faa \
        > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+cnc+oprot2pep106.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        cncRNADB-Translated_ncRNA.8-100.faa \
        openprot-human-openprot-2_2-refprots+altprots+isoforms-min_2_pep-GRCh38.106-uniprot2022_06_01.8-100.faa \
        > $ucat_db/$mfa 
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+oprot2pep106.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        openprot-human-openprot-2_2-refprots+altprots+isoforms-min_2_pep-GRCh38.106-uniprot2022_06_01.8-100.faa \
        > $ucat_db/$mfa 
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+cnc+oprot2pep.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        cncRNADB-Translated_ncRNA.8-100.faa \
        openprot-human-openprot-2_2-refprots+altprots+isoforms-min_2_pep-uniprot2022_06_01.8-100.faa \
    > $ucat_db/$mfa
fi


mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+cnc+cand+sorf.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        cncRNADB-Translated_ncRNA.8-100.faa \
        db-hUniprot_23candidates_20373.8-100.faa \
        db-hUniprot_20373_sorf_15835.8-100.faa \
    > $ucat_db/$mfa
fi

mfa=mabad+mprot+smp2+naug+ribo+spencer+nc+uniprot+ensembl+cnc+cand+sorf.faa
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $ucat_db/$mfa"
    cat db-smprot_openprot_SEPs_MAbad_*.faa \
        MicroProteinDB*.faa \
        SmProt-SmProt2_*circ*.faa \
        SmProt-SmProt2_*highConf*.faa \
        SmProt-SmProt2_*Known*.faa \
        SmProt-SmProt2_*Literature*.faa \
        SmProt-SmProt2_*MS*.faa \
        SmProt-SmProt2_*Ribo*.faa \
        SmProt-SmProt2_*non-AUG* \
        SPENCER-*.faa \
        nc*.faa \
        uniprot*.faa \
        ensembl-Homo_sapiens.GRCh38.pep.abinitio.8-100.faa \
        ensembl-Homo_sapiens.GRCh38.pep.all.8-100.faa \
        cncRNADB-Translated_ncRNA.8-100.faa \
        db-hUniprot_23candidates_20373.8-100.faa \
        db-hUniprot_20373_sorf_15835.8-100.faa \
    > $ucat_db/$mfa
fi



# Finally, make the composite database that we will use as reference
# of all _known_ proteins in the genome, joining openprot, uniprot, 
# ncbi and ensembl
#
mfa=openprot+uniprot+ncbi+ensembl.fasta
if [ ! -s $ucat_db/$mfa ] ; then
    echo "Making $mfa"
    echo \
    cat human-openprot-2_2-refprots+altprots+isoforms-uniprot2022_06_01.faa \
        uniprot-uniprotkb*.faa \
        ncbi-refseq-human*.faa \
        ensembl*pep.all.faa \
        | sed -e '/^>/ s/$/DB=${mfa%.fss}/g' \
        > $ucat_db/$mfa
fi
if [ ! -s ${mfa#*.}_+100.faa ] ; then
    ../bin/select_prot.sh $mfa  # generates openprot+uniprot+ncbi+ensembl_+100.faa
fi


cd ..

exit

#--------------------- OBSOLETE ---------------------
# make non-redundant composite databases
cd $fa_uprot_multi_db
mkdir -p ucat.1


# This will be done in a separate step running remove_redundance.sh
cd $ucat_db

for i in *.faa ; do
    name=${i%.faa}
    ext=${i##*.}
#echo $name ;continue
    if [ ! -s $ucat_db_1/${name}.1.faa ] ; then
#echo $name ;continue
        cd-hit -i $i -o $ucat_db_1/${name}.1.$ext -c 1.0 -M 256000 -T 24 \
               -s 1.0 -S 0 -n 5 -l 8 -t 1 -d 0 -g 1

    fi
done

cd ..


