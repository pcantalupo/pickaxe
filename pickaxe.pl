#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use GetTokens;
use pa;

my $params = {};
# removed node, advres, queue, indexdir=s
GetOptions ($params, "aidfile|a=s",
                      "config=s", "annotconf=s",  # pickaxe config file and annotater config file
                      "walltime|w=s", "cpus|c=i", "email=s",
                      "jobsdir=s",
                      "virusindex=s", "virusfasta=s",   # reference bowtie2 index to align unmapped reads against
                      "subx=s@",     # prepend subtraction index to default set of subx indices
                      "rm_default_index",   # remove default set of subtraction indices
                      "trim_left=i", "trim_right=i",
                      "assembler=s",  # clc_assembler or megahit
                      "skip_repeatmasker",  # skip repeatmasker step
                      "remotetax",    # only get taxonomy information remotely from NCBI (annotater parameter)
                      "collate=s",    # value is 'blast', 'virus' or 'both'
                      "extendcontigs",
                      "stats",        # TBD
                      "exitsubx",   # exit after subx (used when you just want the unmapped reads)
                      "exitassembly", # exit after assembly (used when you just want unmapped and assembled contigs)
                      "exitrm",     # exit after repeatmasker (used when you want unmapped and repeatmasked contigs)
                      "exitent",    # exit after entropy (used when you want the quality contigs: assembly.RM.fa.good.fa.500bp.ent.fa)
                      "negreflist=s",
                      "adapterset=s", "adapterfile=s",   # adapterset not implemented yet in preprocessing.sh
                      "shell|s",
                      "test|t", "force|f",
                      "debug", "help|h",
                      ) or usage();
usage() if (exists $params->{help});
usage() if (exists $params->{aidfile} && ! -e $params->{aidfile});
usage() if (exists $params->{config} && ! -e $params->{config});
usage() if (exists $params->{annotconf} && ! -e $params->{annotconf});

# Set some defaults
$params->{skip_repeatmasker} //= 0;
$params->{remotetax} //= 0;
$params->{exitsubx} //= 0;
$params->{exitassembly} //= 0;
$params->{exitrm} //= 0;
$params->{exitent} //= 0;

# Determine AIDs for processing
$params->{aidfile} = "" if (!exists $params->{aidfile});
my @aid = GetTokens($params->{aidfile});     # each token is an analysis id
usage() if (!@aid);

# Create Pickaxe object
my $pa = pa->new(aids => \@aid, params => $params);

# Run Pickaxe or subcommand
if ($params->{collate}) {
  $pa->collate;
}
elsif ($params->{extendcontigs}) {
  $pa->extendcontigs;
}
elsif ($params->{stats}) {
  $pa->stats;
}
else {
  $pa->runjobs();
}


sub usage {
  print <<USAGE;
PICKAXE version 1.0

Pickaxe identifies potential infectious agents in NGS sequencing reads

Usage:
    pickaxe -a SAMPLEFILE | SAMPLE

 For collating results across samples
    pickaxe.pl --collate both -a SAMPLEFILE | SAMPLE

 For extending contigs for each sample
    pickaxe.pl --extendcontigs -a SAMPLEFILE | SAMPLE

Pickaxe has the ability to search for viruses in various databases and
files.  SAMPLE represents an identifier such as an SRA SRR accession number,
BAM file or FASTQ file.  These identifiers can be put in a SAMPLEFILE (one
identifier per row).

SAMPLE format
The format for SAMPLE is: [ SRR [...] | BAM [...] | FASTQ [...] ]
SRR [...]      NCBI SRA SRR accessions (or ERR or DRR)
BAM [...]      BAM files such as sample1.bam, sample2.bam
FASTQ [...]    FASTQ files.  If you want to combine a set of fastq files
into one group, make sure that each fastq file in the set has the same
prefix and use that prefix as the FASTQ value (example: given foo1.fq and
foo2.fq, then use the prefix foo as the FASTQ value).


Running PICKAXE

First, you need to run pickaxe on an SRR, BAM or FASTQ file(s).  Pickaxe
will put its output in a folder with the same name as the SRR value.  If the
input is a BAM or FASTQ, the output folder name will have contain a
.pickaxe suffix. The output files are as follows.

  aaslurm.job
      the job that was submitted to SLURM or run in the terminal if the
      --shell option was provided
  aaslurm.job.*.slurm.out
      the output from the job
  assembly.RM.fa
      the assembled contigs (raw)
  assembly.ra.tsv
      the relative abundance of the contigs (numreads \t seqLen \t relabund)
  unmapped.fq.gz
      the unmapped reads after performing read subtraction against bowtie2
      indexes such as hg38 
  kv.vrs.hg19.bam
      the alignment results of unmapped reads to the reference database 
  preprocess.log
      the log file of cutadapt and prinseq
  annotator/
      The annotater results directory.  In this directory you will find a
      *BE.report.txt file which contains the BLAST pipeline results and
      annotation for each quality contig (see below for description of
      quality contigs)
  software_versions.tsv
      the versions of software used


Collate results

After running pickaxe, step 2 is to collate the annotater results and
reference database alignment results (from one sample or many) by running
the subcommand 'collate'.  It will generate two files called
report_BLAST.tsv and report_ViralRefSeq.tsv files.

Extending contigs

Sometimes, contigs are part of a longer genome but were not joined during
the assembly process.  Optionally, you can attempt to connect these contigs
by running the subcommand 'extendcontigs'.  The output folder is
'extendedcontigs'.  This processs aligns the viral contigs against the raw
assembled sequences found in assembly.RM.fa using BLASTN for each SAMPLE. 
It performs two rounds of extensions.  For the first extension, ext1.fa can
be empty which means there was no contig extension for that contig.  If a
contig was extended, a second round of BLAST occurs with the ext1 sequence. 
In the second extension, ext2.fa can be empty which means there was no 2nd
contig extension.  If a contig was extended a second time, no further
extension is attempted.


Options:

  --config CONFIGFILE
               Pickaxe configuration file (default ./pickaxe.config)
  --annotconf ANNOTATERCONFIGFILE
               Annotater configuration file (default ./annotater.config)
  --virusindex VINDEX
               Pickaxe requires a Bowtie2 index to which it aligns unmapped
               reads for detection of known viruses in your dataset.  Must
               be an absolute path unless index is found in the path given
               in the env variable BOWTIE2_INDEXES (default                
  --virusfasta VIRUSFASTA
               Absolute path to the fasta file used to generate the VINDEX
  --subx INDEX [--subx ...]
               Specify one or more Bowtie2 indexes to subtract your reads
               against.  Must be an absolute path unless index is found in
               the path given in the environment variable BOWTIE2_INDEXES
  --rm_default_index
               Force Pickaxe to skip subtraction against a set of default
               Bowtie2 indexes
  -w TIME, --walltime TIME
               TIME format is HH:MM:SS (default 48:00:00)
  -c CPUS, --cpus CPUS
               Number of CPU cores to use (default 6)
  --assembler ASSEMBLER
               Supported assemblers are clc_assembler and megahit (default:
               clc_assembler)

Dependendies - listed version may not be absolutely required.  Other
versions were not tested

  samtools	1.2
  bowtie2	2.3.3
  cutadapt	1.18
  PRINSEQ-lite	0.20.4
  RepeatMasker	open-4.0.7
  clc_assembler	5.1.1.184548-201811011136 (or megahit 1.2.9)
  blast+	2.7.1+

Optional

  rapsearch	v2.24
  sra-toolkit	2.9.2



USAGE
  exit 1;
}

