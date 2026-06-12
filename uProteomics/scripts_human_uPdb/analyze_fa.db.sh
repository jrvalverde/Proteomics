#!/bin/bash

# MAKE REFERENCE DB
#   we use openprot-all, but as it is based on uniprot2022, we will
#   concatenate it with uniprot-2026, ncbi-2026 and ensembl-2026,
#   though in principle we only should need uniprot, but just in 
#   case... redundancies, after all, will be removed in the uniqueness
#   clustering step
cd db
mfa=openprot+uniprot+ncbi+ensembl.fasta
if [ ! -s $mfa ] ; then
    echo "Making $mfa"
    cat human-openprot-2_2-refprots+altprots+isoforms-uniprot2022_06_01.faa \
        uniprotkb*.faa \
        ncbi-refseq_human*.faa \
        ensembl*pep.all.faa \
        | sed -e '/^>/ s/$/DB=${mfa%.fss}/g' > $mfa
fi
if [ ! -s ${mfa#*.}_+100.faa ] ; then
    ../bin/select_prot.sh $mfa  # generates openprot+uniprot+ncbi+ensembl_+100.faa
fi
grep -c '^>' *.faa > ../db_SIZES
cd ..

#---------------------------------------------
# PROCESS ALL DATABASES (INCLUDING REFERENCE)
#---------------------------------------------

#   we start by processing all crude databases to get a feeling for their
# actual size, so we may then use this approach to make educated guesses 
# about which databases we can or cannot include
#   in so doing we also process the full proteome (openprot+updates
# and openprot+updates > 100 aa)
#

# make unique full-protein databases
#   we use the full fasta file from cd-hit/usearch because
#   we want only one sequence per cluster to digest, but
#   we want to consider all protein sequences, not only proteins 
#   in size-one clusters
#
echo "Making non-redundant full protein databases"
mkdir -p db.1
bin/remove_redundance.sh db
cd db.1
grep -c '^>' *.faa > ../db.1_SIZES
cd ..

# there is only one unique sequence per cluster
if [ ! -e db_unique ] ; then ln -s db.1 db_unique ; fi
cd db_unique
grep -c '^>' *.faa > ../db_unique_SIZES
cd -

# digest and select peptides between 8 and 30 a.a.
#       SINCE we will only consider peptides 8-30a.a.
# in the microprotein analysis, we can ignore larger
# sizes here as well
#
mkdir -p db_unique_digested
bin/digest.sh db_unique
cd db_unique_digested
grep -c '^>' *.faa > ../db_unique_digested_SIZES
cd ..

if [ ! -e db_pep_8-30 ] ; then ln -s db_unique_digested db_pep_8-30 ; fi
cd db_pep_8-30
grep -c '^>' *.faa > ../db_pep_8-30_SIZES
cd ..

# remove redundancy in peptides 
mkdir -p db_pep_8-30.1
bin/remove_redundance.sh db_pep_8-30
cd db_pep_8-30.1
grep -c '^>' *.faa > ../db_pep_8-30.1_SIZES
cd ..

# extract specific peptides (peptides that appear only in one protein)
# peptide cluster size = 1
mkdir -p db_pep_8-30.1.sp
bin/extract_specific_peptides.03.sh db_pep_8-30.1
cd db_pep_8-30.1.sp
grep -c '^>' *.faa > ../db_pep_8-30.1.sp_SIZES
cd ..

# extract sequences in db.1 that contain specific peptides
# This may have a potential problem:
#   if a protein is encoded more than once in the genome with
# exactly the same sequence, we have considered it only once in
# db_unique, so if it had a specific peptide it would be selected,
# however, if the two copies are indistinguishable, it is as well
# to consider them as one protein for detection purposes
#
#   NOTE THIS WORKS WITH CD-HIT (which creates .clstr files) 
#   FOR USEARCH IT MAKES A PSEUDO .clstr FILE (with unique-only seqs)
mkdir -p db_unique.sp
bin/extract_specific_sequences_from_original_fasta.04.sh db.1 db_pep_8-30.1
cd db_unique.sp
grep -c '^>' *.faa > ../db_unique.sp_SIZES
cd ..


#
# Now repeat allowing for partial digestion
#
#   activate CPU-based AI environment
source ~/contrib/venv/ai-cpu/bin/activate
#   where DeepDIA is installed
INST=~/contrib/DeepDIA

mkdir -p db_ppep_8-30

# NOTE
#   we should TEST pre-processing the fasta file to switch ',' to '%',
#   spaces in the '>' line by '·', and ';' by '!' to see if we can preserve the
#   full description line, e.g.:
#       cat $fasta | sed -e '/>/s/,/%/g' -e '/>/ /·/g' > $preproc_fasta
#   this way, we could recover the full info after DeepDIA processing
#   (DeepDIA apparently preserves only up to the first space) by doing
#   the reverse operation after converting the .obs.50.csv file to
#   fasta.
#

# make digestion files
for f in db.1/*.faa ; do
    fasta=${f}
    name=${fasta%.*}
    name=${name##*/}
    dig_csv=${name}.dig.csv
    if [ ! -s $fasta ] ; then echo "$fasta is empty or missing" ; continue ; fi
    echo "DeepDIA digesting $fasta"
    
    if [ ! -s db_ppep_8-30/$name.pre.faa ] ; then
        # if the preprocessed fasta doesn't exist, then create it
        cat $fasta \
        | sed -e '/>/s/,/%/g' -e '/>/s/ /·/g' -e '/>/s/;/!/g' \
        > db_ppep_8-30/$name.pre.faa
    fi
    # analyze the preprocessed fasta file instead
    fasta=db_ppep_8-30/$name.pre.faa

    if [ ! -s db_ppep_8-30/$name.pep.csv ] ; then
        python $INST/src/digest_proteins.py \
	        --in $fasta \
	        --out db_ppep_8-30/$dig_csv \
	        --group_duplicated \
	        --min_peptide_length 7 \
	        --max_peptide_length 30
    fi
done

# detectable peptides
mkdir -p db_ppep_8-30.detect

for i in db_ppep_8-30/*.dig.csv ; do
    dig_csv=$i
    name=${dig_csv%.*}
    name=${name##*/}
    pred_csv=db_ppep_8-30.detect/$name.obs.csv
    filt_csv=db_ppep_8-30.detect/$name.obs.50.csv
    if [ ! -s $pred_csv ] ; then
        python $INST/src/predict_detectability.py \
	        --in $dig_csv \
	        --model $INST/data/models/detectability/epoch_004.hdf5 \
	        --out $pred_csv
    fi
    if [ ! -s $filt_csv ] ; then
        python $INST/src/filter_peptide_by_detectability.py \
            --peptide $dig_csv \
            --detect $pred_csv \
            --out $filt_csv \
            --min_detectability 0.5
    fi
done
