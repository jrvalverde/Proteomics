db=${1:-db.1}

for i in $db/*.fasta ; do
    fa=${i##*/}
    cl=${db}/${fa}.clstr
    echo $fa $cl
    if [ -s /${db}/${fa} ] ; then continue ; fi
    perl ./extract_Nseq_clstrs.pl $i $cl tmp
    # This will generate one file per fasta sequence in tmp
    # named with a number (0 1 2 3 ... 100 101 ... 1234 ...)
    # we need to concatenate them now
    cd tmp
    mkdir -p ../${db}.1.sp
    ou=../${db}.1.sp/${fa}  # the output fasta file
    rm -f $ou               # remove output if it exists
    truncate -s 0 $ou       # create an empty ouput fasta file
    /bin/ls -v | while read name ; do 
        #cho $name
        cat $name >> $ou    # append to output fasta file
        rm $name
    done
    cd -
    rm ./tmp/*
    echo $i processed
done
