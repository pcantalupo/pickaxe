#!/usr/bin/env perl
use strict;
use warnings;
use PipelineSubx;
use Getopt::Long;
use Data::Dumper;
use File::Copy;

my %o;
my @keys = qw(prefix|p=s threads|t=i seqfile|f=s@ delete_orig_seq|d indxfile|i=s
              keep_unmapped|k keep_subx_sam|s subx|x=s@ rm_default_indx|r help|h
              nosensitive|v nolocal|c
              );
GetOptions (\%o, @keys);

die "Please supply a fastq file (-f|--seqfile)\n" if (! exists $o{seqfile});
#die "Since you have removed the default subx indexes (-r), please supply one or more subx indexes (-x INDX ...)\n" if (! exists $o{subx} && $o{rm_default_indx});

if ($o{help}) {
  print "Usage:\n@keys\n\n";
  exit;
}

if (! exists $o{subx} && $o{rm_default_indx}) {
  my $fastqfile = $o{seqfile}->[0];
  my $totalL = `wc -l $fastqfile | cut -f 1 -d ' '`;
  my $fastqreads = $totalL/4;
  print "No Bowtie2 index supplied\n";
  print "Subtraction pipeline results:\n";
  print "Total reads: $fastqreads\n";
  print "Mapped reads: 0 (0.0)\n";
  print "Unmapped reads: $fastqreads (100.0)\n";

  move($fastqfile, "unmapped.fq");
}
else {
  my $subx = PipelineSubx->new(\%o);
  $subx->run;
}

