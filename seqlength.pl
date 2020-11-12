#!/usr/bin/env perl

# this script outputs the length of each fasta sequence

use strict;
use warnings;
use Bio::SeqIO;

my $file = shift;
my $format = shift || 'fasta';

if (!$file) {
	print "Usage: $0 <file> <format>\n",
			"\tformat can be fasta, qual, ... (see others at http://bioperl.org/wiki/HOWTO:SeqIO\n\n";
	exit;
}

my $in = Bio::SeqIO->new (-file => $file, -format => $format);

while (my $seq = $in->next_seq) {
   
   print $seq->primary_id;
   
   if ($seq->desc) {
      print " " . $seq->desc;
   }
   print "\t", $seq->length, "\n"; 
   
}

