package pa;
use strict;
use warnings;
use File::Path qw(remove_tree);
use File::Basename;
use File::Copy;
use Annotator::Report;
use Annotator::Blast;
use Cwd;
use Carp;
use Data::Dumper;

# Determine PICKAXE directory
my $pickaxe_path = `echo \`which pickaxe.pl\``;
chomp $pickaxe_path;
(my $jobsdir = $pickaxe_path) =~ s/pickaxe.pl$//;
our $JOBSDIR = $jobsdir;

our $PICKAXECONF = "pickaxe.config";
our $ANNOTCONF = "annotater.config";
our $ASSEMBLER = "clc_assembler";
our $OUTDIRSUFFIX = ".pickaxe";
our $VIRUSBAM     = "kv.vrs.hg19.bam";

# Parameters
#    aids - array or hash reference to list of AIDs
# Return - pa object
#
sub new {
  my ($class, %args) = @_;
  my $self = \%args;
  bless $self, $class;
  
  # Load config file
  $self->_load_config;

  if ($self->{params}{'s'}) {
    $self->{params}{shell} = 1;
  }
  $self->{program} = $self->{params}{program};
  $self->{params}{email}    //= "";

  $self->{test}   = $self->{params}{test}   || 0;
  $self->{force}  = $self->{params}{force}  || 0;
  
  $self->{wd}     = getcwd;

  
  
  # Jobs directory
  $self->{params}{jobsdir} //= $JOBSDIR;  
  if ( ! -e $self->jobsdir ) {
    croak "Jobs directory [\"" . $self->jobsdir . "\"] does not exist. Please supply one with '--jobsdir DIR' option\n";
  }

  if ($self->{program} eq 'aaslurm') {
    $self->{params}{collatefile} //= "report_BLAST.tsv";
    $self->{params}{summarizefile} //= "report_ViralRefSeq.tsv";
    $self->{params}{okfile}      //= "aa.OK";
    $self->{params}{contigs}     //= "assembly.RM.fa";
    $self->{params}{jobfile}       = $self->{program} . ".job";
    $self->{params}{trim_left}   //= 0;
    $self->{params}{trim_right}  //= 0;

    # Annotater configuration file
    unless (exists $self->{params}{extendcontigs}) {   # skip if running extendcontigs
      $self->{params}{annotconf} //= $ANNOTCONF;
      if ( ! -e $self->{params}{annotconf} ) {
        croak "Annotater configuration file '" . $self->{params}{annotconf} . "' does not exist. Please supply one with '--annotconf ANNOTCONF' option\n";
      }      
      
      # Set annotconf to absolute path to the annotater configuration file so that annotater (Reann.pl) finds correct file since annotater will be executed in a different directory from working directory
      $self->{params}{annotconf} = Cwd::realpath($self->{params}{annotconf});
      #print "Cwd realpath for 'annotconf': ", Cwd::realpath($self->{params}{annotconf}), $/;
      $self->{params}{annotconffile} = basename $self->{params}{annotconf};
    }

    # virusindex
    #if ( $self->{params}{virusindex} eq "") {
    #  croak "Virusindex not provided.\nEither use option '--virusindex BOWTIE2INDEX' or specify in pickaxe.config configuration file\n\n";
    #}

    if (exists $self->{params}{extendcontigs} || exists $self->{params}{collate}) {
      print "Not checking for valid Bowtie2 index for virusindex since running extendedcontigs or collate\n" if ($self->debug);
    }
    else {
      print "Running aaslurm.job so need to check for Bowtie2 index for virusindex\n" if ($self->debug);

      my $bt2dir = "";
      if (exists $ENV{'BOWTIE2_INDEXES'}) {  $bt2dir = $ENV{'BOWTIE2_INDEXES'} . '/';  }
      $self->{params}{virusindex} //= $bt2dir . "viral.1.1.genomic";
      my $vi = $self->{params}{virusindex};
      
      
      print "Before: virusindex is $vi\n" if ($self->debug);
      if ($vi =~ /^\//) {    # the virus index is a full path already
        print "   no action since full path already\n" if ($self->debug);
      }
      else {
        $vi = $self->{wd} . "/" . $vi;
        print "   adding working directory path to virusindex\n" if ($self->debug);
      
      }
      print "After: virusindex is $vi\n" if ($self->debug);

      $self->{params}{virusindex} = $vi;
      
      my $indxfile  = $vi . ".1.bt2";
      print "indxfile is $indxfile\n" if ($self->debug);
      my $indxfileL = $vi . ".1.bt2l";  # large bowtie2 index
      if ( ! -e $indxfile && ! -e $indxfileL ) {
        croak "Bowtie2 index '$vi' was not found\nEither use option '--virusindex BOWTIE2INDEX' or specify in pickaxe.config configuration file\n\n";
      }
    }
    
    # assembler
    $self->{params}{assembler} //= $ASSEMBLER;
    my $assem = $self->{params}{assembler};
    if ($assem ne "clc_assembler" && $assem ne "megahit") {
      croak "Assembler '$assem' not supported. Please choose either clc_assembler or megahit\n";
    }

    # adapter file
    if (exists $self->{params}{adapterfile} && -e $self->{params}{adapterfile}) {
      my $af = $self->{params}{adapterfile};
      open (my $AF, "<", $af) or croak ("Can't open $af for reading: $!\n");
      while (<$AF>) {
        chomp;
        $self->{params}{adapteroptions} = $_;
      }
      close ($AF);
    }
  }
  else {
    croak "\n$self->{program} is unsupported\n";
  }

  $self->aids;

  print "\nPickaxe configuration\n", Dumper($self),$/ if ($self->debug);

  return $self;
}


sub AUTOLOAD {
  my $self  = shift;
  my $value = shift;

  our $AUTOLOAD;

  ( my $method = lc $AUTOLOAD ) =~ s/.*:://;
  
  if ($value) {
    $self->{params}{$method} = $value;
  }
  return $self->{params}{$method};
}


sub collate() {
  my ($self) = @_;

  ###############
  # TO BE DONE
  # $self->{params}{virusfasta} //= $bt2dir . "viral.1.1.genomic.fna";
  # virus refseqs
  if ( ! exists $self->{params}{virusfasta} || ! -e $self->{params}{virusfasta} ) {
    my $vf = $self->{params}{virusfasta};
    croak "Virusfasta '$vf' was not found\nPlease specify a virus fasta sequence file with '--virusfasta VIRUSFASTA' or specify in pickaxe.config configuration file\n\n";
  }
  ##############


  my $outfile = $self->{params}{collatefile};
  if ( !-e $outfile || (-e $outfile && $self->{force}) ) {
    unlink ($outfile) or die "collate: can't remove $outfile: $!\n" if (-e $outfile);

    #
    # AA program - only collate the hits from individual report files that pass entropy cutoffs
    #
    if ($self->{program} =~ /^aa/ ) {    
      my $collate_method = $self->{params}{collate};
      if ($collate_method eq 'both') {
        $self->_collate_blast();
        $self->_collate_virus();
      }
      elsif ($collate_method eq 'blast') {
        $self->_collate_blast();
      }
      elsif ($collate_method eq 'virus') {
        $self->_collate_virus();
      }
      else {
        croak ("Unknown collate method $collate_method\n");
      }
    }
    else {
      my $g = $self->{refseqs};

      foreach my $aid (sort keys %{$self->aids}) {
        my $bampath = $aid . "/" . $self->{mapdir} . "/" . $self->{finalbam};
        my $command = "summarizebam_by_ref.pl -a $aid -g $g -f $bampath >> $outfile";
        `$command`;
      }
    }
  }
}

sub _collate_blast {
  my ($self) = @_;
  print "Running collate blast\n" if ($self->debug);
  
  my $header = "aid\tseqID\tseq\tseqLength\tpid\tcoverage\te\taccession\tdesc\ttype\tfamily\tspecies\tgenome\talgorithm\tdb\tqstart\tqend\tsstart\tssend\tnsf\tqent\tqshp_ent\tshsp_ent\tshsp_%lc\tra\n";
  my $outfile = $self->{params}{collatefile};
  open (my $out, ">", $outfile) or die "Can't open $outfile: $!\n";
  print $out $header;

  my $rpthash;
  foreach my $aid (sort keys %{$self->aids}) {
    print "\t$aid\n" if ($self->debug);

    my $aidtype = $self->_aid_type($aid);
    if ($aidtype eq 'bam' || $aidtype eq 'fastq') {
      $aid .= $OUTDIRSUFFIX;
    }
    my $reportfile = $aid . "/annotator/ann.wTax.BE.report.txt";

    if (-e $reportfile) {
      my $ar = Annotator::Report->new(report => $reportfile);
      my $passedentr = $ar->pass_filters(use_report => 1);

      my $blastfilter = 0;
      my $ab;
      if (`grep 'qc 90' $aid/annotator/$self->{params}{annotconffile}`) {
        $ab = Annotator::Blast->new( blast => [ "$aid/annotator/ann.0.0.blastn","$aid/annotator/ann.0.1.blastn" ] );
        $ab->runfilter( qc => 50, pid => 80, evalue => 1e-5 );
        $blastfilter = 1;
      }

      foreach (@$passedentr) {
        chomp;
        next if (/seqID/);
        my ($contig, @rest) = split(/\t/, $_);

        if ($blastfilter && $rest[7] eq '') {                  # rest[7] is the 'type' field
          next if ($ab->passfilter($contig));                  # skip over those unassigned sequences that would pass the blast filter criteria
        }

        push (@rest, "NULL");
        $rpthash->{$aid}{$contig} = \@rest;
      }
    }

    my $relabundance = $aid . "/assembly.ra.tsv";
    if (-e $relabundance) {
      open (my $RA, "<", $relabundance);
      while(<$RA>) {
        chomp;
        my ($contig, undef, undef, $ra) = split(/\t/, $_);
        if (exists $rpthash->{$aid}{$contig}) {
          $rpthash->{$aid}{$contig}[-1] = $ra;
        }
      }
    }

    foreach my $contig (keys %{$rpthash->{$aid}}) {
      (my $aid_prefix = $aid) =~ s/\.pickaxe//;
      print $out join("\t", $aid_prefix, $contig, @{$rpthash->{$aid}{$contig} }), "\n";
    }
  }
}

sub _collate_virus {
  my ($self) = @_;

#  my $g = $self->virusfasta;
#  if (! -e $g) {
#    croak "Virus fasta sequence file, $g, not found\n";
#  }

  print "Running collate virus\n" if ($self->debug);

  my $outfile = $self->{params}{summarizefile};
  my $header = "aid\taccession\tvname\tvlength\tseqcov\tavgdepth\taligns\tavgMAPQ\tavgScore\tavgEditDist\n";
  open (my $out, ">", $outfile) or die "Can't open $outfile: $!\n";
  print $out $header;
  close ($out);

  my $virusfasta = $self->virusfasta;
  foreach my $aid (sort keys %{$self->aids}) {
    print "\t$aid\n" if ($self->debug);
    
    my $aidtype = $self->_aid_type($aid);
    if ($aidtype eq 'bam' || $aidtype eq 'fastq') {
      $aid .= $OUTDIRSUFFIX;
    }

    my $bampath = $aid . "/" . $VIRUSBAM;
    $aid =~ s/\.pickaxe//;
    my $command = "summarizebam_by_ref.pl -a $aid -g $virusfasta -f $bampath >> $outfile";
    `$command`;
  }
}


sub extendcontigs {
  my ($self) = @_;
  
  print "Running extend contigs\n" if ($self->debug);

  my $blastreportfile = $self->{params}{collatefile};
  if (! -e $blastreportfile) {  #   "report_BLAST.tsv"
    croak "Report BLAST file '$blastreportfile' not found...exiting\n";
  }
  
  foreach my $aid (sort keys %{$self->aids}) {
    print "\t$aid\n" if ($self->debug);
    
    my $aidtype = $self->_aid_type($aid);
    if ($aidtype eq 'bam' || $aidtype eq 'fastq') {
      $aid .= $OUTDIRSUFFIX;
    }

    chdir($aid);

    my $dir = "extendedcontigs";
    mkdir ($dir) unless (-d $dir);

    chdir ($dir);

    my $contigs = $self->{params}{contigs};    
    (my $id = $aid) =~ s/\.pickaxe//;
    my $command = "run_extendcontigs.sh $id ../../$blastreportfile ../$contigs > extend.out";
    `$command`;
    print $command,$/ if ($self->debug);

    chdir($self->{wd});
  }
}

sub stats {
  my ($self) = @_;

  print "Calculating stats for pickaxe runs\nTO BE FINISHED\n";
}

# Get/Set AIDs
#
# Stores list of AIDs as a hash reference for efficiency
#
# Parameters:
#    aids - can be a hash or array reference that contains AIDs
# Return: hash ref to the AIDs
sub aids {
  my ($self, $aids) = @_;

  my $aidfile = "aid.txt";

  # Different cases for setting / refactoring aids:
  # 1. 'aids' was set in client as an array reference ( aids => \@aid )
  # 2. A reference to aids (array or hash) were passed into this method. Copy over existing aids.
  # 3. Client created PA object with empty AID array reference or
  #    supplied an AID reference to this method. In either event, we do not have a valid HASH reference
  #    of AIDS so try to get AIDs from default aid.txt file

  # case #1
  if (ref $self->{aids} eq "ARRAY" && @{ $self->{aids} } ) { # convert existing array of aids to a hash ref
    $self->{aids} = {  map { $_ => 1 } @{ $self->{aids} }  };
  }
  # case #2
  elsif ($aids) {
    my %tmp;
    if (ref $aids eq "ARRAY") {
      %tmp = map { $_ => 1 } @$aids;     # convert array ref to hash
    } elsif (ref $aids eq "HASH") {
      %tmp = %$aids;
    }
    if (keys %tmp) {
      $self->{aids} = \%tmp;
    }
  }

  # case #3
  if (ref $self->{aids} ne "HASH") {
    $self->{aids} = _get_aids_from_file($aidfile);
  }

  return $self->{aids};
}


sub runjobs {
  my ($self) = @_;
  my $jobfile = $self->{params}{jobfile};
  my $okfile = $self->{params}{okfile};

  foreach my $aid (sort keys %{$self->aids}) {
    my $dir = $aid;
    my $aid_type = $self->_aid_type($aid);
    if ($aid_type eq 'bam' || $aid_type eq 'fastq') {
      $dir .= $OUTDIRSUFFIX;
    }
    mkdir ($dir) unless (-d $dir);
    chdir ($dir);
  
    print $aid;

    # check if OK file exists, if not there is work to be done!
    if (-e $okfile && !$self->{force}) {
      print "\tno work to be done...skipping";
    }
    else {
      my $wd = $self->{wd};
      
      # Negative reference list
      if ( (   $aid_type eq 'aid' || $aid_type eq 'srr' || $aid_type eq 'bam' )
            && exists $self->{params}{negreflist}
         ) {
        my $nrf = $self->{params}{negreflist};
        `ln -s $wd/$nrf`;
      }

      # Link to BAM or FASTQ files in parent directory since the jobfile expects these files in the current directory
      if ($aid_type eq 'bam' || $aid_type eq 'fastq') {
        # if AID represents a single BAM or FASTQ file, then link to it
        if (-e "$wd/$aid") {
          `ln -s $wd/$aid` unless (-l $aid);
        }
        
        # else AID represents a FASTQ prefix, so need to link to each FASTQ file in the group
        else {
          my @fastq_files = glob("$wd/$aid*");
          foreach my $fq (@fastq_files) {
            next if ($fq =~ /$OUTDIRSUFFIX$/);
            `ln -s $fq` unless (-l basename($fq));
          }
        }        
      }

      my $path = $self->jobsdir . "/$jobfile";
      open (my $rawjob, "<", $path) or die "Can't open jobfile $path for reading: $!\n";
      open (my $newjob, ">", $jobfile) or die "Can't open jobfile $jobfile for writing: $!\n";

      my $jobname = "$jobfile." . $aid;         # aaslurm.job.SRR1234      
      my $joboutputfile = "$jobname.out";       # aaslurm.job.SRR1234.out
      if ( -e $joboutputfile ) {
        my $backupfile = $joboutputfile . ".1";       # Careful: this only backups the last output file
        move($joboutputfile, $backupfile);
      }
      
      $self->{params}{walltime} = '48:00:00' if (!exists $self->{params}{walltime});
      $self->{params}{cpus}     = 6          if (!exists $self->{params}{cpus});
      $self->{params}{adapterset} = ""       if (!exists $self->{params}{adapterset});
      $self->{params}{adapteroptions} = ""   if (!exists $self->{params}{adapteroptions});

      while (<$rawjob>) {
        my $email = $self->email; s/__EMAIL__/$email/;
        s/__WALLTIME__/$self->{params}{walltime}/;
        s/__CPUS__/$self->{params}{cpus}/;
        s/__MEM__/$self->{params}{cpus}*16/e;
        s/__OUTPUTFILE__/$joboutputfile/;
        s/__NAME__/$jobname/;
        s/__OKFILE__/$okfile/;
        s/__AID__/$aid/;
        s/__AIDTYPE__/$aid_type/;
        s/__ADAPTERSET__/$self->{params}{adapterset}/;
        s/__ADAPTEROPTIONS__/$self->{params}{adapteroptions}/;
        s/__EXITSUBX__/$self->{params}{exitsubx}/;
        s/__EXITASSEMBLY__/$self->{params}{exitassembly}/;
        s/__EXITRM__/$self->{params}{exitrm}/;
        s/__EXITENT__/$self->{params}{exitent}/;
        s/__TRIMLEFT__/$self->{params}{trim_left}/;
        s/__TRIMRIGHT__/$self->{params}{trim_right}/;
        s/__ASSEMBLER__/$self->{params}{assembler}/;
        
        my $vi = $self->virusindex; s/__VIRUSINDEX__/$vi/;

        if (exists $self->{params}{annotconf}) {
          $ANNOTCONF = $self->{params}{annotconf};
        }
        s/__ANNOTATERCONF__/$ANNOTCONF/;

        if (exists $self->{params}{negreflist}) {
          my $nrl = $self->{params}{negreflist};
          s/__NEGREFLIST__/$nrl/;
        }
        else {
          s/__NEGREFLIST__//;
        }

        if (/__SUBTRACTION__/) {
          my $xparam = "";
          if (exists $self->{params}{subx}) {
            $xparam = join (" -x ", " ", @{ $self->{params}{subx} } );
          }
          if (exists $self->{params}{rm_default_index}) {
            $xparam .= " -r";
          }
          s/__SUBTRACTION__/$xparam/;
        }
        
        print $newjob $_;
      }
      close ($newjob) or die "Can't close jobfile $jobfile: $!\n";

      my $jobid = "00000";
      unless ($self->{test}) {
        unlink $okfile;
        
        my $command = "bash";
        if ($self->{params}{shell}) {
          $command = 'bash';
        }
        elsif ($self->{program} =~ /slurm/i) {
          $command = 'sbatch';
        }

        if ($command ne 'bash') {
          chomp ($jobid = `$command $jobfile`);
        }
        else {
          `bash $jobfile > $joboutputfile 2>&1`;
        }
      }
      print "\t$jobid";
    }

    chdir ($self->{wd});
    print "\n";
  } # done with this AID
}


sub _get_aids_from_file {
  my ($file) = @_;

  my %aids;
  if (-s $file) {
    open (my $fh, "<", $file) or die "Found $file file but cannot open it: $!\n";
    %aids = map { chomp; $_ => 1 } <$fh>;
  }

  return \%aids;
}

sub _load_config {
  my ($self) = @_;

  $self->{params}{config} //= $PICKAXECONF;
  my $cf = $self->{params}{config};
  if (-e $cf) {
    print "Processing configuration file $cf\n" if ($self->debug);
    open (my $fh, "<", $cf);
    while (<$fh>) {
      chomp;
      next if (/^#/);
      my ($name, $value) = split (/=/, $_);    # format: NAME=VAL
      if ($name && $value) {
        $self->$name($value) unless $self->$name;  # don't set the value if it is already set
      }
    }
    close ($fh);  
  }
}

sub _aid_type {
  my ($self, $aid) = @_;

  return 'bam' if ($aid =~ /\.bam$/);
  return 'aid' if ($aid =~ /\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}/);
  return 'srr' if ($aid =~ /^.RR\d+$/);
  return 'syn' if ($aid =~ /^syn\d+$/);
  return 'fastq';
  #return 'unsupported';
}




1;
