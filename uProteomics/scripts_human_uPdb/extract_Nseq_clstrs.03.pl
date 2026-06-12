#!/usr/bin/perl

#note you have to use "-d 0" in the cd-hit run
#note it's better to use "-g 1" in the cd-hit run

# this script reads the .clstr file, it generates a seperate fasta file
# for each cluster of certain fixed size
# for example, if you run cd-hit -i db -o dbout -c 0.6 -n 4 -d 0 -g 1
# then you will have a dbout.clstr
# now run this script:
#    ./make_multi_seq.pl db dbout.clstr multi-seq 20
# to get all clusters of size 20
# if no size is specified, then 1 is assumed (i.e. only print
# fasta sequences for clusters containing only one sequence)
# 

my $fasta = shift;
my $clstr = shift;
my $output = shift;
my $size_cutoff = shift;

die unless (-e $fasta);
print "FASTA $fasta\n";
die unless (-e $clstr);
print "CLUSTER $clstr\n";
die unless ($output);
print "OUTPUT $output\n";
$size_cutoff = 1 unless ($size_cutoff);

if (not -e $output) {
    if ($size_cutoff == 1) {my $cmd = `touch $output`;}
    else                   {my $cmd = `mkdir $output`;}
}

# open cluster file
open(TMP, $clstr) || die;
my %id2cid=();
my $cid = "";
my @ids = ();
my $no = 0;

$k = 0;
# process cluster file line by line
while($ll = <TMP>) {
  # if last line starts with >, then we are at a new cluster entry
  if ($ll =~ /^>/) {
    # before processing this new cluster entry, we check how many
    # members did te previous cluster have ($no) and if it matches
    # the cutoff, we add all the members to the previous cluster
    #if ($no >= $size_cutoff) {
    if ($no == $size_cutoff) {
    #if ($no > 0) {
      print "$k ($no)\n"; $k++;
      # add all cluster ids to the id2cid conversion table and
      # associate them with their cid
      foreach $i (@ids) { $id2cid{$i} = $cid; print "$i: $id2cid{$i} = $cid\n";}
    }

    # once we have saved the previous cluster, we can proceed
    # to get the newly found one
    if ($ll =~ /^>Cluster (\d+)/) { $cid = $1; }
    else { die "Wrong format $ll"; }
    # since we are in a new cluster, there are no members yet
    # so we clear the member ids and set the count to 0
    @ids = ();
    $no = 0;
  }
  # if line does not start with '>' then we are looking at a
  # cluster member of the last detected cluster
  else {
    # append it to an array of current cluster members
    if ($ll =~ /(aa|nt), >(.+)\.\.\./) { push(@ids, $2); $no++; }
    else { die "Wrong format $ll"; }
  }
}
close(TMP);
    # add last members of the cluster
    if ($no == $size_cutoff) {
      foreach $i (@ids) { $id2cid{$i} = $cid;}
    }


foreach $i (@ids) { print "$id2cid{$i}\n" ; }
#print 'DONE';
#exit;

# now that we have a list of seqs associated with their cluster name
# we can go over the fasta file and save each sequence to its
# corresponding cluster fasta file
#
# open the file with the fasta sequences
open(FASTA, $fasta) || die;
# is the output file open?
my $outfile_open = 0;
my $save_seq = 0;

# if the cutoff is 1, we will only save one file: we can open it here
if ($size_cutoff == 1) {
  open(OUT, ">> $output") || die "can not open file to write $output/$cid";
  $outfile_open = 1;
}

# for each line in the fasta file
while($ll=<FASTA>) {
  # if it is a new entry
  if ($ll=~ /^>(\S+)/) {
    # get the ID
    my $id = $1;
    # lookup the $id in the id2cid table
    my $cid = $id2cid{$id};
    # if found, we need to save the header and all subsequent lines until
    # a new entry is found
    if (defined($cid)) {
      if ($size_cutoff != 1) {
        # ensure we save it to the correct fasta file
        close(OUT) if ($outfile_open);
        open(OUT, ">> $output/$cid") || die "can not open file to write $output/$cid";
        $outfile_open = 1;
      }
      $save_seq = 1;
    }
    else {
      # if not found in the id-to-cluster table, do not save this sequence
      $save_seq = 0;
    }
  }
  if ($save_seq) {print OUT $ll;}
}
close(FASTA);
close(OUT) if ($outfile_open);


