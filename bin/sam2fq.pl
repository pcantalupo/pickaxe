#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $negref_file = "";
my $flag_keep   = -1;
my $format      = "fastq";
my $paired = 0;
my $prefix = "foo";         # foo_1.fq, foo_2.fq
GetOptions("n=s"  => \$negref_file,
            "f=i" => \$flag_keep,
            "t=s" => \$format,
            "paired" => \$paired,
            "prefix=s" => \$prefix,
            );

# process negative reference sequence name list
my %nr = ();
if ($negref_file) {
   open(IN, "<", $negref_file) or die "sam2fq: can't open file $negref_file: $!\n";
   while (<IN>) {
      chomp;
      $nr{$_}++;   
   }
}


my %preads;   # cache paired reads until you have both reads of a pair so you can output them in order
my ($one, $two, $error);
if ($paired) {
   open ($one, ">", $prefix . "_1.fq");
   open ($two, ">", $prefix . "_2.fq");
   open ($error, ">", $prefix . "_error.fq");
}

while (<>) {
   next if /^@/;
   chomp;
   my ($id, $flag, $rname, $seq, $qual) = (split /\t/, $_)[0,1,2,9,10];

   my $print = 0;
   if ($flag_keep == -1) {   # accept all flag values
      $print = 1 unless (exists $nr{$rname});    # print only if refname doesn't appear in negative refname list
   } elsif ($flag_keep & $flag) {
      $print = 1;
   } elsif (%nr) {
      $print = 1 unless (exists $nr{$rname});    # keep seqs that are not mapped to the negative refname list
   }

   next unless ($print == 1);    # skip this sequence unless you want to print it

   if ($flag & 16) {
      $seq = reverse($seq);
      $seq =~ tr/[ACGTacgt]/[TGCAtgca]/;
      $qual = reverse($qual);
   }

   if ($paired) {
      my $toPrint = ($format eq "fastq") ? join("","@",$id,"\n",$seq,"\n","+","\n",$qual,"\n") : join("",">",$id,"\n",$seq,"\n");
      if ($flag & 64) {
         $preads{$id}{1} = $toPrint;
      }
      elsif ($flag & 128) {
         $preads{$id}{2} = $toPrint;
      }
      else {
         print $error $toPrint;
      }
      
      if ($preads{$id}{1} && $preads{$id}{2}) {
         print $one $preads{$id}{1};
         print $two $preads{$id}{2};
         delete $preads{$id};
      }
   }
   else {      
      if ($format eq "fastq") {
         print "@", $id,   "\n",
               $seq,       "\n",
               "+",        "\n",
               $qual,     "\n";
      }
      else {
         print ">", $id, "\n",
               $seq,     "\n";
      }
   }
}

# for all the paired reads that didn't have a mate in the BAM file, need to
# print the remaining reads to Error file (i.e.  singletons)
if ($paired) {
   foreach my $id (keys %preads) {
      foreach (keys %{$preads{$id}}) {
         print $error $preads{$id}{$_};
      }
   }
}

