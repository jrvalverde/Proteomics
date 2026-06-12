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

# PROCESS ALL DATABASES (INCLUDING REFERENCE)
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
#   NOTE THIS WORKS WITH CD-HIT (which creates .clstr files) BUT NOT USEARCH
mkdir -p db_unique.sp
bin/extract_specific_sequences_from_original_fasta.04.sh db.1 db_pep_8-30.1
cd db_unique.sp
grep -c '^>' *.faa > ../db_unique.sp_SIZES
cd ..



# * * * * * * * * * * * * * * * * * * * * 
#   PROCESS MICROPROTEIN DATABASES
# * * * * * * * * * * * * * * * * * * * * 


# make microprotein databases directory
echo "Selecting microproteins (size <= 100)"
mkdir -p uprots
bin/select_microprot.sh
cd uprots
echo "Calculating sizes"
grep -c '>' * > ../uprots_SIZES
cd ..

# make non-redundant microprotein databases
# if two microprotein entries have exactly the same
# sequence we will consider them as the same microprotein
#   -- this means we need to use the .fasta files of either
#   cd-hit or u-clust, not just the unique (cluster size=1) ones
echo "Making non-redundant microprotein databases"
mkdir -p uprots.1
bin/remove_redundance.sh uprots
cd uprots.1
echo "Calculating sizes"
grep -c '^>' *.faa > ../uprots.1_SIZES
cd ..

echo "FULLY digest non-redundant microprotein databases"
mkdir -p uprots.1_digested
bin/digest.sh uprots.1 uprots.1_digested
cd uprots.1_digested
echo "Calculating sizes"
grep -c '^>' *.faa > ../uprots.1_digested_SIZES
cd ..
ln -s uprots.1_digested uprots.1_pep_8-30
cp uprots.1_digested_SIZES uprots.1_pep_8-30_SIZES

echo "Making augmented microprotein+protein FULLY digested peptide databases"
mkdir -p uprots.1_pep_8-30+
for i in uprots.1_digested/*.faa ; do
    name=`basename $i`
    name=${name%.faa}
    echo "Preparing $i"
    if [ ! -s uprots.1_pep_8-30+/${name}.faa ] ; then
        cat $i \
            db_pep_8-30/openprot+uniprot+ncbi+ensembl_+100.faa \
            > uprots.1_pep_8-30+/${name}.faa
    fi
done
cd uprots.1_pep_8-30+
echo "Calculating sizes"
grep -c '^>' *.faa > ../uprots.1_pep_8-30+_SIZES
cd ..


echo "Remove redundancy in FULLY digested microprotein peptides"
mkdir -p uprots.1_pep_8-30+.1
bin/remove_redundance_usearch.sh uprots.1_pep_8-30+
cd uprots.1_pep_8-30+.1
grep -c '^>' *.faa > ../uprots.1_pep_8-30+.1_SIZES
cd ..

echo "Extract microproteins that have unique fully digested peptides"
mkdir -p uprots.1.sp
bin/extract_specific_sequences_from_original_fasta.05.sh \
        uprots.1 \
        uprots.1_pep_8-30+.1
cd uprots.1.sp
grep -c '^>' *.faa > ../uprots.1.sp_SIZES
cd ..



# * * * * * * * * * * * * * * * * * * * * 
#       COMBINE MICROPROTEIN DATABASES
# * * * * * * * * * * * * * * * * * * * * 

# make concatenated databases from non-redundant microprotein databases
# these may have many copies of the same microprotein because it may
# have been present in different databases with different entry names,
echo "Building concatenated databases"
bin/make_composite_dbs.sh
cd ucat
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat_SIZES
cd ..

# since the same protein may have been present in different databases
# with different entry names we cluster and select all cluster representatives
# (all .fasta files in cd-hit and usearch)
mkdir ucat.1
bin/remove_redundance.sh ucat
cd ucat.1
grep -c '^>' *.faa > ../ucat.1_SIZES
cd ..

# digest the databases to obtain peptide databases
echo "Digesting concatenated databases"·
mkdir -p ucat.1digested
bin/digest.sh ucat.1 ucat.1_digested
# create convenience soft link
if [ ! -L ucat.1_pep_8-30 ] ; then
    ln -s ucat.1_digested ucat.1_pep_8-30
fi
cd ucat.1_pep_8-30
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30_SIZES
cd ..

# make non-redundant peptide databases
echo "Making non-redundant peptide databases"
mkdir -p ucat.1_pep_8-30.1
bin/remove_redundance.sh ucat.1_pep_8-30
cd ucat.1_pep_8-30.1
echo "Calculating sizes (all peptides)"
grep -c '^>' *.faa > ../ucat.1_pep_8-30.1_SIZES
cd ..

if [ "NOT" = "YES" ] ; then
# Filter peptides by their detectability using a Random Forest
# trained on a hundred+ parameters. 
# This is VEEEEEEEERY SLOOOOOOOOW, so we will only do it on
# a few chosen small databases
mkdir -p ucat.1_pep_8-30.1_scored
mkdir -p ucat.1_pep_8-30.1_detectable
bin/pep_detectability.sh ucat.1_pep_8-30.1
cd ucat.1_pep_8-30.1_detectable
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30_detectable.1_SIZES
cd ..

# once finished, we can link back the filtered databases
# into ucat.1_pep_8-30.1 to add them to the downstream
# processing steps
for i in ucat.1_pep_8-30.1_detectable/* ; do
    name=`basename $i`
    name=${name%.faa}
    cp $i ucat.1_pep_8-30.1/${name}_obs.faa
done

fi


# Now that we have the microprotein peptide databases, we
# can filter for uniqueness by combining them with the full
# PROTEIN proteome (WITHOUT microproteins) peptide 
# database: any common peptides between both or within any 
# will be deemed non-specific and considered a target for 
# removal.
#
# Note that we use the database '+100.faa' which contains proteins > 100a.a.
mkdir -p ucat.1_pep_8-30+
for i in ucat.1_pep_8-30/*.faa ; do
    name=`basename $i`
    name=${name%.faa}
    echo "Preparing $i"
    if [ ! -s ucat.1_pep_8-30+/${name}+.faa ] ; then
        cat $i \
            db_pep_8-30/openprot+uniprot+ncbi+ensembl_+100.faa \
            db_pep_8-30/openprot+uniprot+ncbi+ensembl_+100.faa \
            > ucat.1_pep_8-30+/${name}.faa
    fi
done
cd ucat.1_pep_8-30+
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30+_SIZES
cd ..

# cluster the databases to get redundant and non-redundant peptides
#   Any peptide in the microprotein database that is also in the
#   protein > 100aa database, or in the protein > 100aa database that
#   is in also in the microprotein database will cluster in a size > 1 
#   cluster and we will get size=1 clusters with all microprotein and 
#   protein unique peptides
#
#   We get only microprotein unique peptides by adding twice the 
#   protein > 100aa database, because then all protein peptides would 
#   be duplicated and size=1 clusters would only be feasible from micro-
#   proteins, but this would duplicate cd-hit running time, and it is
#   not strictly necessary if we later extract only microproteins with 
#   unique peptides
#
mkdir -p ucat.1_pep_8-30+.1
bin/remove_redundance.sh ucat.1_pep_8-30+
cd ucat.1_pep_8-30+.1
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30.1+.1_SIZES
cd ..


# Make databases with only the sequences that contain specific unique
#   peptides: since we know peptides are unique, if a peptide belongs
#   to a protein > 100 a.a. then its full sequence will not be in the
#   microprotein database and since we only extract sequences from the
#   microprotein databases, only microproteins with unique microprotein
#   peptides will be selected, and not any protein > 100aa
#
#   We use a specific modified version of the extraction script that
#   takes into account the special names of the database+reference
#   files
#
echo "Creating uniquely identifiable full sequence databases"
mkdir -p ucat.1.sp
# using the full peptide databases
# CHECK:  NOTE THIS WORKS WITH CD-HIT (which creates .clstr files) BUT NOT USEARCH
bin/extract_specific_sequences_from_original_fasta.05.sh ucat.1 ucat.1_pep_8-30+.1
cd ucat.1.sp
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1.sp_SIZES
cd ..



# * * * * * * * * * * * * * * * * * * * * 
#   REMOVE ADDITIONAL SEQUENCES
# * * * * * * * * * * * * * * * * * * * * 


# To substract specific sequences we first concatenate them
# to the selected ucat.1.sp databases TWICE and then we remove all
# duplicates. The reason for adding them twice is that if we do
# add them only once, then common seqs would be removed, but
# non-common seqs would be added, by having them twice we make
# sure that all sequences in the 'substract' sets will be repeated
# and hence removed.
make ucat.1.sp+substract
for i in ucat.1.sp/mabad* ; do
    db=${i#*/}
    echo $db
    name=${db%.faa}
    echo $name
    cat $i \
        substract/*.faa \
        substract/*.faa \
    > ucat.1.sp+substract/${name}+substract.faa
done

bin/remove_redundance.sh ucat.1.sp+substract

