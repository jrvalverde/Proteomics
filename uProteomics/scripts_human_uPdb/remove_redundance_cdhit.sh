db=${1}

cd $db

mkdir -p ../${db}.1

for i in ./*.faa ; do
    echo "Deduplicating $i"
    if [ -s ../${db}.1/$i ] ; then continue ; fi
    date
    # -i input sequence
    # -o output
    # -c sequence identity threshold
    # -M memory limit
    # -T number of threads
    # -s length difference cutoff (0..1)
    # -S length difference cutoff in amino acid from cluster representative
    # -n word length (default 5)
    # -l length of throwaway sequences !!!
    # -d length of description in .clstr file (0 stop at firsdt space)
    # -g use global identity (1) or local identity (0)
    time cd-hit -i $i \
        -o ../${db}.1/${i} \
        -c 1.0 \
        -M 256000 \
        -T 24 \
        -s 1.0 \
        -S 0 \
        -n 5 \
        -l 8 \
        -t 1 \
        -d 0 \
        -g 1
    date
    # the fasta file contains all cluster representatives
    # independent of cluster size
    # to extract the representative of clusters with size=1
    # we need to use the .clstr files to extract them
done |& tee ../${db}.1/uniq-cd-hit.log
cd -

exit
# extract_unique_sequences will look for fasta files and
# generate .faa files with representatives of size=1 
# clusters
#   that script expects all-cluster files to be named .fasta
#bin/extract_unique_sequences.03.sh ${db}
#

exit


db=${1}

mkdir -p ${db}.1
cd $db

for i in ./*.faa ; do
    echo "Deduplicating $i"
    if [ -s ../${db}.1/$i ] ; then continue ; fi
    date
    time cd-hit -i $i -o ../${db}.1/$i -c 1.0 -M 256000 -T 24 \
          -s 1.0 -S 0 -n 5 -l 8 -t 1 -d 0 -g 1
    date
done |& tee ../${db}.1/uniq-cd-hit.log
cd -
