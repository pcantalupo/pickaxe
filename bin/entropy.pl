#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use SeqUtils;
use Getopt::Long;

sub usage {

  my $out = <<HEREDOC;

Usage: $0 [-e] -f FILE

Script requires a fasta formatted file as input. The default action is to output a fasta file containing sequences that have entropy value greater than the min entropy value (> 65). If you supply the -e option, the script will output a row of three columns for each sequence: entropy value   \t   fasta id    \t   sequence

-m     Minimum entropy value (default: 65)
-e     Entropy only 


HEREDOC

  print $out;
  exit 0;
}

my $minE = 65;
my $eonly = 0;
my $file = '';
my $help = 0;
GetOptions ("min|m=i"  => \$minE,
            "eonly|e"  => \$eonly,
            "file|f=s" => \$file,
            "help|h"   => \$help,
            ) or exit;
&usage if ($help);
die "Please supply fasta file (file|f)\n" if (!$file || ($file ne '-' && !-e $file));

my $seqio = ($file eq '-') ? Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta') : 
                             Bio::SeqIO->new(-file => $file);

while (my $seq = $seqio->next_seq) {

  my $ent = entropy ($seq->seq);
  if ($eonly) {
    print $ent,"\t",$seq->primary_id,"\t",$seq->seq,"\n";
  }
  else {
    if ($ent > $minE) {
      print ">" . $seq->primary_id . "\n" . $seq->seq . "\n";
    }
  }
}



