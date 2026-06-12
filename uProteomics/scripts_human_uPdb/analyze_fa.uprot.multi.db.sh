
cd ..
topdir=`pwd`
cd -

full_pep_dbs=$topdir/fa.db/db_pep_8-30


# * * * * * * * * * * * * * * * * * * * * 
#       COMBINE MICROPROTEIN DATABASES
# * * * * * * * * * * * * * * * * * * * * 

# make concatenated databases from non-redundant microprotein databases
# these may have many copies of the same microprotein because it may
# have been present in different databases with different entry names,
echo "Building concatenated databases"
cd ..
bin/make_composite_dbs.sh
cd -
cd ucat
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat_SIZES
cd ..

# since the same protein may have been present in different databases
# with different entry names we cluster and select all cluster representatives
# (all .fasta files in cd-hit and usearch)
mkdir ucat.1
../bin/remove_redundance.sh ucat
cd ucat.1
grep -c '^>' *.faa > ../ucat.1_SIZES
cd ..


#------------------------------------
#   SIMULATE FULL DIGESTION WITH R
#------------------------------------

# digest the databases to obtain peptide databases
echo "Fully digesting concatenated databases"·
mkdir -p ucat.1_digested
../bin/digest.sh ucat.1 ucat.1_digested
# create convenience soft link
if [ ! -L ucat.1_pep_8-30 ] ; then
    ln -s ucat.1_digested ucat.1_pep_8-30
fi
cd ucat.1_pep_8-30
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30_SIZES
cd ..

# make non-redundant peptide databases
echo "Making non-redundant fully digested peptide databases"
mkdir -p ucat.1_pep_8-30.1
bin/remove_redundance.sh ucat.1_pep_8-30
cd ucat.1_pep_8-30.1
echo "Calculating sizes (all peptides)"
grep -c '^>' *.faa > ../ucat.1_pep_8-30.1_SIZES
cd ..

if [ "$RF_DETECT" = "TRUE" ] ; then
    # Filter peptides by their detectability using a Random Forest
    # trained on a hundred+ parameters. 
    # This is VEEEEEEEERY SLOOOOOOOOW, so we will only do it on
    # a few chosen small databases (mabad+*)
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
            $full_pep_dbs/openprot+uniprot+ncbi+ensembl_+100.faa \
            > ucat.1_pep_8-30+/${name}+.faa
    fi
done
cd ucat.1_pep_8-30+
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30+_SIZES
cd ..

# cluster the databases to get redundant and non-redundant peptides
#   Any peptide in the microprotein database that is also in the
#   protein > 100aa database, or in the protein > 100aa database that
#   is in also in the microprottein database will cluster in a size > 1 
#   cluster and we will get size=1 clusters with all microprotein and 
#   protein unique peptides
#
#   We could get only microprotein unique peptides by adding twice the 
#   protein > 100aa database, because then all protein peptides would 
#   be duplicated and size=1 clusters would only be feasible from micro-
#   proteins, but this would duplicate cd-hit running time, and it is
#   not necessary if we later extract only microproteins with unique
#   peptides
#
mkdir -p ucat.1_pep_8-30+.1
bin/remove_redundance.sh ucat.1_pep_8-30+
cd ucat.1_pep_8-30+.1
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1_pep_8-30+.1_SIZES
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
echo "Creating fully digested uniquely identifiable microprotein databases"
mkdir -p ucat.1.sp
# using the full peptide databases
# NOTE THIS WORKS WITH CD-HIT (which creates .clstr files) BUT
# WITH USEARCH WE DO NOT GET CLUSTER FILES, SO WE MAKE A *FAKE* ONE
bin/extract_specific_sequences_from_original_fasta.05.sh ucat.1 ucat.1_pep_8-30+.1
cd ucat.1.sp
echo "Calculating sizes"
grep -c '^>' *.faa > ../ucat.1.sp_SIZES
cd ..


#===========================================================================
#   WORK OUT PARTIALLY DIGESTED DETECTABLE PEPTIDE DATABASES WITH DEEPDIA
#===========================================================================
#
#
# We start from ucat.1
#
# digestion peptides (including partial digestion) will go into ucat.1_ppep_8-30
#
#   activate CPU-based AI environment
source ~/contrib/venv/ai-cpu/bin/activate
#   where DeepDIA is installed
INST=~/contrib/DeepDIA

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
mkdir -p ucat.1_ppep_8-30
for f in ucat.1/*.faa ; do
    fasta=${f}
    name=${fasta%.*}
    name=${name##*/}
    dig_csv=${name}.dig.csv
    if [ ! -s $fasta ] ; then echo "$fasta is empty or missing" ; continue ; fi
    echo "DeepDIA digesting $fasta"
    
    if [ ! -s ucat.1_ppep_8-30/$name.pre.faa ] ; then
        # if the preprocessed fasta doesn't exist, then create it
        cat $fasta \
        | sed -e '/>/s/,/%/g' -e '/>/s/ /·/g' -e '/>/s/;/!/g' \
        > ucat.1_ppep_8-30/$name.pre.faa
    fi
    # analyze the preprocessed fasta file instead
    fasta=ucat.1_ppep_8-30/$name.pre.faa

    if [ ! -s ucat.1_ppep_8-30/$name.pep.csv ] ; then
        python $INST/src/digest_proteins.py \
	        --in $fasta \
	        --out ucat.1_ppep_8-30/$dig_csv \
	        --group_duplicated \
	        --min_peptide_length 7 \
	        --max_peptide_length 30
    fi
done

# detectable peptides
mkdir -p ucat.1_ppep_8-30.detect

for i in ucat.1_ppep_8-30/*.dig.csv ; do
    dig_csv=$i
    name=${dig_csv%.*}
    name=${name##*/}
    pred_csv=ucat.1_ppep_8-30.detect/$name.obs.csv
    filt_csv=ucat.1_ppep_8-30.detect/$name.obs.50.csv
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
    # do final postprocessing from obs.50.csv when converting to fasta file
    # reverting the preprocessing step
done
