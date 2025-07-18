#!/usr/bin/env perl

# Process RepeatMasker masked fasta sequences in two steps:
# 1. Remove sequences that do not contain a stretch of 50 consecutive base pairs (i.e. GATC)
# 2. Remove sequences that contain more than or equal to 40% N's.

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $stretch = 50;   # must have at least a run of this many nucleotides to pass the 'Stretch' criteria
my $percent = 40;   # must have less than this percent of Ns in the sequence to pass the 'Percent Ns' criteria

my $fasta;
my $mask = 'N';       # default is to look for N's in the sequence (can set to 'tcga' to look for lowercase tcga's)
GetOptions ("fasta|f=s" => \$fasta,
            "mask|m=s"  => \$mask,
            );
die "No fasta provided" if !$fasta;
my $name = $fasta;
$name =~ s/\.[^\.]*$//;

my $seqI = new Bio::SeqIO(-file => $fasta, -format => 'fasta');

# output filehandle for the quality sequences that pass the criteria as a quality sequence
my $seqO = new Bio::SeqIO(-fh => \*STDOUT, -format => 'fasta');

my $num_pass_stretch = 0;
my $num_quality      = 0;
my $total_seqs       = 0;
while(my $seq = $seqI->next_seq){
	$total_seqs++;
	if($seq->seq =~ /[^$mask]{50,}/){
		$num_pass_stretch++;
	
		my @mm = split '', $seq->seq;
		my $Ncount = 0;
		foreach (@mm){
			if($_ =~ /[$mask]/){
				$Ncount++;
			}
		}
		
		if($Ncount/$seq->length < .4) {
			$num_quality++;
			$seqO->write_seq($seq);   # only outputs seqs that passed both criteria (Quality seqs)
		}
	}
}

my $num_nopass_stretch = $total_seqs - $num_pass_stretch;
my $num_nopass_percent = $num_pass_stretch - $num_quality;
my $percent_quality    = $num_quality/$total_seqs * 100;

print STDERR
"Num seqs >=",$stretch," nucleotide stretch: $num_pass_stretch\n",
"Num seqs <",$stretch," nucleotide stretch : $num_nopass_stretch\n",
"\n",
"Num seqs >=",$percent,"% ",$mask,"'s              : $num_nopass_percent\n",
"\n",
"Total Quality Seqs: $num_quality (", $percent_quality, "%)\n",
"Total Input Seqs  : $total_seqs\n";
		

