# * * * * * * * * * * * * * * * * * * * * 
#   PROCESS MICROPROTEIN DATABASES
# * * * * * * * * * * * * * * * * * * * * 


# make microprotein databases directory
echo "Selecting microproteins (size <= 100)"
mkdir -p uprots
cd ..                   # must be run from the top level
bin/select_microprot.sh
cd -                    # return to where we were
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
