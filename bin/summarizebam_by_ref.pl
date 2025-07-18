#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;

my $bamfile = '';
my $genome  = '';
my @annots  = ();        # these are additional annotations that will get prepended to output
                         # example: -a LIHC,93b7  --> the first two columns of output will always be LIHC\t93b7
my $pass = 0;
my $minuniqpos = 150;
my $minalign   = 50;
my $ignoredups = 0;
my $remove = 0;          # remove reads that have the same alignment score in the 2nd best hit (XS field)
my @tags = qw/AS NM XS/;    # tracking AlignmentScore and NumberMismatches
GetOptions('annots|a=s{,}' => \@annots,
          'file|f=s'       => \$bamfile,
          'genome|g=s'     => \$genome,
          'ignoredups|d'   => \$ignoredups,
          'minuniqpos|u=i' => \$minuniqpos,
          'minalign|n=i'   => \$minalign,
          'pass|p'         => \$pass,
          'remove|r'       => \$remove,
          );
if ($bamfile eq ""){ print "Please supply '-f BAMFILE' on command line\n"; exit; }
if ($genome eq "") { print "Please supply '-g GENOME' on command line\n"; exit; }


#
# process genome file
#
my %ref = ();
my $seqio = Bio::SeqIO->new(-file => $genome);
while (my $seq = $seqio->next_seq) {
  my $rid = $seq->primary_id;
  $ref{$rid}{name} = $seq->description;
  $ref{$rid}{len}  = $seq->length;
}


#
# process BAM file
#
# get number of alignments, align score and some quality information
baminfo(\%ref, $bamfile, $remove);

# get samtools depth information for file
refdepth(\%ref, $bamfile);

#print Dumper(\%ref),"\n";
#exit;

#
# output
#
foreach my $rid (keys %ref) {
  if (exists $ref{$rid}{nalign}) {
    my $rlen = $ref{$rid}{len};
    
    my $avgdepth = $ref{$rid}{td} / $rlen;          # average depth per base for genome    
    my $upos = $ref{$rid}{tp};
    my $seqcov   = $upos / $rlen * 100;    # percent of genome length that is covered by at least one read

    my $nalign  = $ref{$rid}{nalign};
    my $avgmapq = $ref{$rid}{mapq}/$nalign;

    my @tag_avg = ();
    foreach my $t (@tags) {
      next if ($t eq 'XS');
      push(@tag_avg, $ref{$rid}{tags}{$t}{total} / $ref{$rid}{tags}{$t}{count});
    }
        
    print join("\t", @annots, $rid, $ref{$rid}{name}, $rlen, $seqcov, $avgdepth, $nalign, $avgmapq, @tag_avg);
    if ($pass) {
      my $toPrint = ($nalign >= $minalign && $upos >= $minuniqpos) ? "\t1" : "\t0";
      print $toPrint;
    }
    print "\n";
  }
}



#
# SUBROUTINES
#
sub refdepth {
  my ($ref, $file) = @_;

#  my @depth = `samtools depth $file`;
  open (my $in, "samtools depth $file | ");
  while (<$in>) {    # each row in depth output is a successive nucleotide position in the reference sequence
    chomp;
    
    my ($rid, $s, $d) = split(/\t/, $_);
    $$ref{$rid}{tp}++;           # total positions
    $$ref{$rid}{td} += $d;       # total depth
  }
}


sub baminfo {
  my ($ref, $file, $remove) = @_;
#  my @samlines = `samtools view $file`;

  open(my $in, "samtools view $file | ");
  while (<$in>) {
    chomp;
    my ($qname, $flag, $rid, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @rest) = split (/\t/,$_);
    next if ($rid eq '*');    # unmapped read
    next if ($ignoredups && $flag & 1024);

    my $keep = 1;
    my %tagvalue;
    foreach (@rest) {
      if (/(\w+):i:(.+)$/) {
        my $name = $1;
        my $value = $2;
        foreach (@tags) {
          if ($name eq $_) {
            $tagvalue{$name} = $value;
          }
        }
      }
    }

    if ($remove && exists $tagvalue{'XS'} && $tagvalue{'AS'} == $tagvalue{'XS'}) {
      $keep = 0;
    }

    if ($keep) {
      foreach my $t (keys %tagvalue) {
        $$ref{$rid}{tags}{$t}{count}++;
        $$ref{$rid}{tags}{$t}{total} += $tagvalue{$t};
      }
      $$ref{$rid}{nalign}++;
      $$ref{$rid}{mapq} += $mapq;
    }
  }
}

