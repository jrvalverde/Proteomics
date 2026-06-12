CSV=$1
name=${CSV%.csv}
echo "$CSV -> $name.faa"
if [ -s $name.faa ] ; then echo "--> $name already done, exiting" ; exit ; fi
# skip header line and reverse any preprocessing:
# copy anything up to first comma and anything up to second comma
#   and output '>second copy (description)\nfirst copy (peptide sequence)
# then, in '>'-containing lines change '%' to ',', '·' to ' ' and '!' to ';'
# NB: this last step can also be done with "tr '%·!' ', ;' but since we have
#   already loaded sed, we can just reload it and make the change specific 
#   to '>'-lines
#   note: tr(1) should normally be faster than sed(1), but since the operating
#   system uses lazy discarding of memory, and sed has just been used, we
#   expect the second sed(1) to exploit the just-used, still-in-memory, copy
#   while tr(1) would require a file system search, disk read and memory
#   load
#   note: since we are likely to use this script repeatedly, in modern
#   computers with lots of memory, all programs (tail, sed, tr) will likely
#   be in cache after the first run, so there won't be a difference in
#   speed likely, but sed allows us to ensure changes apply only to '>'-lines
#   providing a safety net
tail -n +2 $1 \
| sed -e 's/^\([^,]\+\),\([^,]\+\),.*/>\2\n\1/g' \
| sed -e  '/>/s/%/,/g' -e '/>/s/·/ /g' -e '/>/s/!/;/g' \
> $name.faa

