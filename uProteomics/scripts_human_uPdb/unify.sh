file=$1
name=${file%.*}
ext=${file##*.}

cd-hit -i $file -o $name.1.$ext -c 1.0 -M 256000 -T 24 \
          -s 1.0 -S 0 -n 5 -l 8 -t 1 -d 0 -g 1

