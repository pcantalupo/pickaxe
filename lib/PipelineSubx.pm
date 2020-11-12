package PipelineSubx;
use strict;
use warnings;
use Carp;
use POSIX;
use File::Basename;

our @indx = qw/NT_rRNA.all Homo_sapiens.GRCh37.74.ncrna Homo_sapiens.GRCh37.74.cdna.all
                       refseq.human.rna              ucsc_genes_hg38
                       hg19                          hg38
                       human_genomic.00              human_genomic.01
                       human_genomic.02              human_genomic.03
                       human_genomic.04              human_genomic.05
                       human_genomic.06              human_genomic.07
                       human_genomic.08              human_genomic.09
                       human_genomic.10              human_genomic.11
                       human_genomic.12
                       /;
sub new {
  my $class = shift;
  my $self = shift;

  $self->{prefix} //= "unmapped";
  $self->{threads} //= 2;
  $self->{bt2_dir} = $ENV{'BOWTIE2_INDEXES'} . '/';
  
  bless $self, $class;
  $self->initialize;  

  return $self;
}

# Initialize - build array to hold the subtraction indexes that will be used in the pipeline
sub initialize {
  my ($self) = @_;

  # three types of subx indexes
  # 1) default set
  if ($self->{rm_default_indx}) {
    @indx = ();
  }
  
  # 2) indexes specified by user on command line (put them to the front of the line)
  if ( $self->{subx} ) {
    unshift (@indx, @{$self->{subx}});
  }
  
  # 3) indexes specified by user in a file (put them to the front of the line)
  if (exists $self->{indxfile}) {
    my $file = $self->{indxfile};
    open (my $in, "<", $file) or die "Can't open indxfile $file in initialize method of PipelineSubx: $!\n";
    my @x;
    while (<$in>) {
      chomp;
      push (@x, $_);
    }
    unshift (@indx, @x);
  }
  
  foreach my $indx (@indx) {
    # check for index (both regular and large index) in local and in BOWTIE2_INDEXES directory
    my $indxfile = $indx . ".1.bt2";
    my $indxfileL = $indx .".1.bt2l";  # large bowtie2 index
    if (   ! -e $indxfile
        && ! -e $self->{bt2_dir} . $indxfile
        && ! -e $indxfileL
        && ! -e $self->{bt2_dir} . $indxfileL) {
      croak "a Bowtie2 index '$indx' was not found (search paths: ./ and $self->{bt2_dir}\n";
    }
  }
  
  if (!@indx) {
    croak "No subx indexes specified\n";
  }

  # Output version info
  print "Bowtie2 index version information (index, numseqs, date, filesize)\n";
  foreach my $indx (@indx) {
    my ($indx_path, $indxfile) = $self->_get_indx_files($indx);
    my $lastseq_info = `bowtie2-inspect -s $indx_path | tail -n1`;
    my ($N) = $lastseq_info =~ /^Sequence-(\d+)/;
    open (my $fh, "<", $indxfile) or die ("Can't open $indxfile for reading: $!\n");
    my @stats = stat($fh);
    my $S = $stats[7];
    my $D = strftime("%m/%d/%Y", localtime($stats[9]));   # index 9 is modified time (see perldoc -f stat)
    print join("\t", $indx_path, $N, $D, $S), $/;
    close ($fh);
  }
}

# Run - run the subtraction pipeline
#         1) map reads to subtraction index
#         2) calculate metrics for mapping
#         3) get unmapped reads from SAM file
sub run {
  my ($self) = @_;
   
  my $i = 0;
  for (my $i = 0; $i <= $#indx; $i++) {
    $self->{i} = $i;
    $self->mapper;
    $self->mapmetrics;
    $self->unmapped;        
  }

  $self->report;
}


# For an index such as `foo`, returns two scalars
# 1. the full path to the index (indx_path)                    : /PATH/TO/foo
# 2. the full path to the '.1.bt2' or '.1.bt2l' file (indxfile): /PATH/TO/foo.1.bt2
sub _get_indx_files {
  my ($self, $indx) = @_;
  
  my $indx_path = "";
  my $indxfile = "";
  my @suffix = ("", "l");
  foreach (@suffix) {
    my $indx_suffix = ".1.bt2" . $_;
    $indxfile = $indx . $indx_suffix;
    if (-e $indxfile) {
      $indx_path = $indx;
      last;
    }
    elsif (-e $self->{bt2_dir} . $indxfile) {
      $indxfile = $self->{bt2_dir} . $indxfile;
      $indx_path = $self->{bt2_dir} . $indx;
      last;
    }
  }

  return ($indx_path, $indxfile);
}


sub mapper {
  my ($self) = @_;
  my $i = $self->{i};
  my $indx = $indx[$i];
  
  my ($indx_path) = $self->_get_indx_files($indx);
  croak "Path to $indx not found\n" if ( ! $indx_path );   # this shouldn't occur since there is a check for this in 'initialize' method
  
  foreach my $seqfile (@{$self->{seqfile}} ) {
    croak "$seqfile does not exist\n" if (! -e $seqfile);
  }
  my $seqfiles = join(",", @{$self->{seqfile}} );
  $indx = basename($indx);
  $self->{sam} = "subx.$indx.sam";
  $self->{logfile} = "subx.$indx.log";
  
  my $sensitive = "--very-sensitive";
  if (exists $self->{nosensitive}) {
    $sensitive = "";
  }
  
  my $local = "--local";
  if (exists $self->{nolocal}) {
    $local = "";
  }
    
  my $cmd = "bowtie2 $sensitive $local -p $self->{threads} -x $indx_path -U $seqfiles -t > $self->{sam} 2> $self->{logfile}";
  print "$self->{i}: $cmd\n";
  `$cmd`;

  if ($i == 0 && $self->{delete_orig_seq}) {   # default is to keep original sequence files (i.e. foo_1.fastq and foo_2.fastq)
    unlink (@{$self->{seqfile}});
  }
   
  if ($i > 0 && ! exists $self->{keep_unmapped}) {  # default is to delete subsequent input sequences
    unlink (@{$self->{seqfile}});
  }
}

# Mapmetrics - parse Bowtie2 log file for the number of reads aligned, num read mapped
#                 and num reads unmapped
#
sub mapmetrics {
  my ($self) = @_;
  my $i = $self->{i};
  my $indx = $indx[$i];

  open (my $in, "<", $self->{logfile});

  my ($total, $unmapped);  
  while (<$in>) {
    chomp;
    if (/(\d+) reads; of these:/) {
      $total = $1;
    }  
    if (/(\d+) \(\S+\) aligned 0 times/) {
      $unmapped = $1;
    }  
  }
  my $m = { indx => $indx, total => $total, unmapped => $unmapped, mapped => ($total - $unmapped) };
  push (@{$self->{metrics}}, $m);
  
  close ($in);
  unlink ($self->{logfile});
}

sub unmapped {
  my ($self) = @_;
  my $i = $self->{i};
  my $indx = $indx[$i];

  # determine unmapped sequence filename
  $indx = basename($indx);     # only want basename if 'indx' contains directory paths
  my $unmapped = ($self->{i} == $#indx) ? ($self->{prefix}) : "unmapped.$indx";
  $unmapped .= ".fq";
 
  # extract unmapped from SAM file
  my $cmd = "sam2fq.pl -f 4 $self->{sam} > $unmapped";
  #print "$i: $cmd\n";
  `$cmd`; 
  
  # delete SAM by default
  unless (exists $self->{keep_subx_sam}) {
    unlink ($self->{sam});
  }
  
  $self->seqfiles( [ $unmapped ] );
}

# Report - output number of reads being aligned, mapped and unmapped for each subtraction step
#        - also, reports the starting total number of reads, how many were mapped and unmapped
#           for entire subtraction pipeline 
sub report {
  my ($self) = @_;
  
  my $i = 0;
  my $T;
  my $M = 0;
  print "Subtraction pipeline results:\n";
  foreach (@{$self->{metrics}}) {
    if ($i == 0) {
      $T = $_->{total};
    }
    $M += $_->{mapped};
    print join("\t", $_->{indx}, $_->{total}, $_->{mapped}, $_->{unmapped}),"\n";
    $i++;
  }
  print "Total reads: $T\n";
  my $Mper = $M/$T*100;
  printf ("Mapped reads: %i (%.1f)\n", $M, $Mper);
  my $U = $T-$M;
  my $Uper = $U/$T*100;
  printf "Unmapped reads: %i (%.1f)\n", $U, $Uper;
}


# seqfiles - get/setter function
#
# seqfile - array reference of one or more sequence files                                  
sub seqfiles {
  my ($self, $seqfiles) = @_;
  
  if ($seqfiles) {
    $self->{seqfile} = $seqfiles;
  }  
  return $self->{seqfile};
}
