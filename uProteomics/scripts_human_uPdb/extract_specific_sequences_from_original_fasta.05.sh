#!/bin/bash
#
#   For generic clusters we need the original fasta, because the
# clustered fasta only contains the representative sequence, but
# not the other cluster sequences, but for specific peptides 
# (peptides that belong only to one protein and therefore form 
# size-1 clusters) we can save time using the clustered fasta because
# we only need the centroid/representative in it.
#
#   We can also save time if we prepare things so that we can run 
# many scripts in parallel in separate machines.
#
#   We can further improve speed by changing extract_Nseq_clstrs.pl
# so that when size_cutoff is 1, the output destination is taken 
# not as a directory to save clusters separately, but as a file to
# save all size-1 clusters together. THIS IS INCONSISTENT BEHAVIOUR
# COUNTERINTUITIVE AND DANGEROUS, SO THE MODIFIED SCRIPT SHOULD BE
# HEAVILY DOCUMENTED, BUT THE SPEEDUP GAINED CAN BE IMPRESSIVE.
#
#   IN THIS VERSION WE DO USE THE MODIFIED extract_Nseq_clstrs.03.pl
# THAT EXTRACTS THE SEQUENCE NAME (WITHOUT _CLID) SO IT MATCHES THE
# ORIGINAL FASTA SEQUENCE NAME
#

# reserve __ (two underscores) for comments and make it read-only
unset __ ; declare -r __ # and now we can use ${__+ comment } or ${__# comment}

db=${1:-db_unique.1}            # the directory with full fasta sequences
clstr_dir=${2:-db_pep_8-30.1}   # the directory with unique clusters

# Next is a specific version which, instead of using the ID output
# to the
extract_Nseq_clstrs=bin/extract_Nseq_clstrs_from_original_fasta.04.pl

# create output directory
mkdir -p ${db}.sp

# filter each fasta in $db using the matching cluster file in $clstr_dir
for i in $db/*.faa ; do
    # get the fasta file name (without directory path)
    fa=${i##*/}             # fasta file
    fan=${fa%.faa}          # fasta file name
    # generate the cluster file name
    cl=${clstr_dir}/${fa}.clstr             # cluster file name
    # define the output file name
    ou=${db}.sp/${fa}                     # the final output fasta file
    echo $fa $cl $ou
    if [ ! -s "$cl" ] ; then
        # we must be using usearch12, which does not generate a cluster
        # file (neither a .uc file due to bugs in the current beta)
        # so, we need to create a mock one: we cannot create a full
        # one because we do not know all members in multi-sequence
        # (size > 1) clusters, only the representative, but since we 
        # only are interested in specific (size=1) clusters, and since
        # we know the requirements of $extract_Nseq_clusters, we can 
        # create a PSEUDO-cluster file that contains
        #   >Cluster #
        #   0   ....aa, >fasta-entry-until-first-space... *
        # for each sequence in the fasta file that has ;size=1;
        #
        # change cluster file name to reflect the origin
        cl=${cl/clstr/uclstr}
        if [ ! -s "$cl" ] ; then
            echo "    creating $cl"
            #echo "
            cat ${clstr_dir}/$fa           ${__#cat file contents } \
            | sed -e 's/;size=/ ;size=/g'  ${__# insert a space before ;size= } \
            | tr '\n' '@'                  ${__# make all content one line } \
            | sed -e 's/@>/\n>/g'          ${__# insert a \n before each sequence entry } \
            | grep ';size=1;'              ${__# select entries with cluster size=1 } \
            | nl -n ln                     ${__# add left-justified line numbers } \
            | sed -e 's/^\([0-9]\+\)\s\+/@Cluster \1\n/g' ${__# insert "Cluster " before line numbers at start of line } \
            | sed -e '/>/s/\([^[:space:]]\+\).*/0\tNNNaa, \1... */g' ${__# insert '0\tNNNaa, ' before entry name } \
            | sed -e 's/@/>/g'             ${__# switch back remaining @ to \n } \
            > $cl
            #"
        fi
    fi
    #continue
    #exit
    if [ -s "$ou" ] ; then continue ; fi  # already done
    
    # ensure the output file exists to avoid race conditions
    # if run in parallel
    truncate -s 0 "$ou"       # create an empty ouput fasta file
    
    echo "    filtering $i"
    echo "        using $cl"
    echo "           to $ou"
    perl $extract_Nseq_clstrs "$i" "$cl" "${ou}" 1
    # NOTES:
    # We exploit here the default cluster size setting (4th argument), 
    # which is 1.
    # This makes use of a specific version of our perl script that
    # distinguishes when size_cutoff is one, and then, instead of
    # using separate files for each cluster, saves all size-1 clusters
    # to a single file
done
