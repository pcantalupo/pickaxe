#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Bamtools;
use Bio::SeqIO;

sub usage {
  print "
$0 -c CONTIGS -r READS [--cpus 4 --cleanup -h -f FORMAT]

Supply the reads (fastq or fasta) that were used with the assembler to
create the CONTIGS file (fasta).

--cpus     number of threads (default: 4)
--cleanup  remove BAM file and Bowtie index that is created

Output from the script:
1) Files (can be deleted with --cleanup option)
  a) bowtie2 index of the CONTIGS file
  b) BAM file of READS mapping to the CONTIGS
2) STDOUT:
  a) table of relative abundances of the contigs. Note: if a contig does not
  have any alignments in the BAM file, it will be absent from this
  relabuncance table.
3) STDERR:
  a) bowtie2 output of alignment statistics

";
exit;
}

my $format = "fastq";
my ($contigs, $reads, $help) = ("","",0);
my $cleanup;
my $cpus = 4;
GetOptions ("c|contigs=s" => \$contigs,
            "r|reads=s"   => \$reads,
            "f|format=s"  => \$format,
            "h|help"      => \$help,
            "cleanup"     => \$cleanup,
            "t|cpus=s"    => \$cpus,
            ) or usage;

usage if ($help);
usage if (!$contigs || !$reads);

#
# create bowtie index
my $base = basename($contigs);
my ($prefix, $suffix) = $base =~ m/(.*)\.(.+)$/;
`bowtie2-build -f $contigs $prefix`;

#
# map reads to assembly
my $bamout = $reads . ".bam";
my $format_option = "";   # default to Fastq option for Bowtie2
if ($format eq 'fasta') {
  $format_option = "-f";
}
my $cmd = "bowtie2 -p $cpus $format_option -U $reads -x $prefix | samtools view -Sb - | samtools sort - $reads";
#print $cmd, $/;
`$cmd`;
my @refaligns = bam2refnumaligns("$bamout");


my %contiglen;
my $in = Bio::SeqIO->new(-file => $contigs,-format => 'fasta');
while (my $seq = $in->next_seq) {
  $contiglen{$seq->primary_id} = $seq->length;
}

my $unmapped;
foreach my $refaligns (@refaligns) {
  my ($aligns, $ref) = split(/\t/, $refaligns);
  next if ($ref eq '*');
  
  my $len = $contiglen{$ref};
  if ($len <= 0) {
    print "$len is <= 0 for $ref. Refaligns line is <$refaligns>\n";
    exit;
  }
  my $relabundance = $aligns/$len; 
  print join("\t", $ref, $aligns, $len, $relabundance), "\n";
  
  delete $contiglen{$ref};
}

# The remaining keys in %contiglen are those contigs that do not appear in
# the BAM file because there were no reads mapping to the contig by bowtie. 
# Somehow the assembler is creating a set of contigs where none of the reads
# map back to the contig.
foreach my $ref (keys %contiglen) {
  print join("\t", $ref, "0", $contiglen{$ref}, "0"), "\n";
}

=head
# Determine the QNAMEs (i.e. reads) that aligned to each Reference (i.e. contigs)
my %refs = get_reads_for_each_qname( $bamout );
open (my $REFS_TO_READS, ">", $prefix . ".refs2qname.tsv");
foreach (sort keys %refs) {
  print $REFS_TO_READS
        $_,
        "\t",
        join( ",", @{$refs{$_}} ),
        $/;
}
close ($REFS_TO_READS);
=cut

if ($cleanup) {
  my @bt2 = glob($prefix . "*bt2");
  unlink (@bt2, $bamout);
}
