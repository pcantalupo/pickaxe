#!/usr/bin/env perl

# This script reads in fasta sequence identifiers and searches a 
# fasta sequence file for those sequences. If found, it outputs the 
# fasta sequence.

use strict;
use warnings;
use GetTokens;
use Bio::SeqIO;
use Getopt::Long;

my ($firstword, $help, $seqfile, $seqidfile);
GetOptions ("firstword|f"    => \$firstword,
            "help|h"         => \$help,
            "seqidfile|i=s"  => \$seqidfile,
            "seqfile|s=s"    => \$seqfile,
            ) or usage();

&usage if (!$seqfile);
&usage if $help;

my %seqids = map {$_ => 1} GetTokens($seqidfile);

my $seqin = Bio::SeqIO->new(-file => $seqfile, -format => 'fasta');
my $seqout = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'fasta');

my %seqs = ();
while (my $seq = $seqin->next_seq() ) {
   # build the fasta header
   my $fasta_header = $seq->display_id;
   if ($seq->desc) {
      $fasta_header .=  " " . $seq->desc;
   }

   my $seqid;
   if ($firstword) {
      ($seqid) = $fasta_header =~ /^(\S+)/;   # capture the first word of the fasta header
   } else {
      $seqid = $fasta_header;    # use entire fasta header as the seqid
   }

   # check to see if this sequence is one that we want and if so, output it
   if (exists $seqids{$seqid}) {
      $seqout->write_seq($seq);
   } 
}


sub usage {
   print <<EOF;

Usage:
$0 [-f] [-i SEQIDfile] -s SEQS_FILE SEQID [ SEQID..]

     -i SEQIDfile  contains a unique identifier on each line
                   to search for in the SEQS_FILE
     -s SEQS_FILE  must contain fasta sequences where the sequence
                   is on one line
     -f            Use first word in Fasta header as the sequence id


     Since the script uses GetTokens for capturing the sequence ids
     from a pipe, you cannot pipe a sequence file into the script. 
     For example, -s /dev/stdin will not work but redirection does
     work (i.e. < seqs.fa) 
     
EOF

   exit;

}
