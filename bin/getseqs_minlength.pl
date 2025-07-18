#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;


my $file = "";
my $minlen = "";
GetOptions ('f|file=s' => \$file,
        'l|len=i' => \$minlen,
        );

if ($file eq "" || $minlen eq "") {
  print "Usage: $0 -f FASTAFILE -l MINLENGTH\n\n";
  exit;
}

my $seqio;
if ($file eq '-') {
  $seqio = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta');
}
else {
  $seqio = Bio::SeqIO->new(-file => $file);
}
my $seqout = Bio::SeqIO->newFh(-fh => \*STDOUT, -format => 'fasta');

my $tl = 0;
my $nseqs = 0;
while (my $seq = $seqio->next_seq ) {

  my $len = $seq->length;
  if ($seq->length >= $minlen) {
    $nseqs++;
    $tl+= $len;
    print $seqout $seq;
  }
}

print STDERR "Total length of $nseqs seqs that are >= $minlen: $tl\n";

