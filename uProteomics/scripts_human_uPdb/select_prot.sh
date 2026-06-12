in=${1}

name=${in%.*}
ext=${in#*.}

if [ "${ext}" = "faa" ] ; then echo "Must NOT be a .faa file" ; exit ; fi

out=${name}_+100.faa

min=101
max=-1

echo "Extracting sequences > 100 a.a. from"
echo "    $in"
echo "into"
echo "    $out"
seqkit seq -g -m $min -M $max "$in" > "$out"
