#!/usr/bin/env perl

package ExtendSequence;
use Carp;
use Bio::SeqIO;

sub new {
  my $class = shift;
  my $self = shift;
  
  bless $self, $class;
  
  if (! -e $self->{query} ) {
    croak "$self->{query} file does not exist\n";
  }
  $self->{queryseq} = Bio::SeqIO->new(-file => $self->{query})->next_seq->seq;   # get just the sequence of the query fasta sequence
  
  $self->{pid} //= "90";
  
  $self->prepare_blast_db;
  
  return $self;
}


sub prepare_blast_db {
  my ($self) = @_;
  
  $self->{blastdb} //= "tmp";
  $self->{target} //= "assembly.RM.fa";
  
  if (! -e $self->{target}) {
    croak "$self->{target} file does not exist\n";
  }
  
  if ( ! -e $self->{blastdb} . ".nhr" && ! -e $self->{blastdb} . ".nal" ) {
    `makeblastdb -dbtype nucl -in $self->{target} -out $self->{blastdb}`;
  }
}

sub run {
  my ($self) = @_;
  $self->{blastout} //= "extendblast.bls";
  my $q = $self->{query};
  my $d = $self->{blastdb};
  `blastn -task blastn -query $q -db $d -outfmt '6 std qlen slen' > $self->{blastout}`;  
}

sub report {
  my ($self) = @_;
  
  open (my $B, "<", $self->{blastout}) or croak ("Can't open $self->{blastout}: $!\n");
  
  while (<$B>) {
    chomp;
    
    $self->{currenthit} = $_;
    my $eseq = $self->_extendseq();
    
    if ($eseq) {
      print $eseq;
    }    
  }
}

sub _extendseq {
  my ($self) = @_;
    
  my ($q, $s, $pid, $hsplen, undef, undef, $qs, $qe, $ss, $se, $e, undef, $qlen, $slen) = split /\t/, $self->{currenthit};
  
  return undef if ($q eq $s || $pid < $self->{pid} || $hsplen < $self->{minoverlap});


=head
Blast hit is at Start of Query - 4 possibilites
                   valid     qs  qe  ss  se
start subject       no       1   15  1   15
start subject RC    yes      1   15  15  1
end   subject       yes      1   15  486 500
end   subject RC    no       1   15  500 486
=cut

  my $eseq;
  my $subjseq = $self->_getsubjseq($s);

  # Start Q
  # 	extended seq will be S->Q
  #
  if ( $qs == 1 && $se eq 1 ) {  # Start S (RC) hits query       TESTED
    $eseq = ">${s}rc-$q overlap " . ($slen - ($ss - $se) ) . " $slen\n";
    $eseq .= substr($self->_getsubjseq($s)->revcom->seq, 0, $slen - ($ss - $se + 1) ) . $self->{queryseq} . $/;
  }
  
  elsif ( $qs == 1 && $se == $slen ) { # End S hits query          TESTED
    if ($se - $ss + 1 == $subjseq->length()) {   # occurs when Subject seq exists completely at start of Query seq
      #print STDERR "Subject $s exists completely at 5' end of $q\n";
      return undef;
    }

    $eseq = ">$s-$q overlap $ss $se\n";
    $eseq .= substr($self->_getsubjseq($s)->seq, 0, $ss - 1) . $self->{queryseq} . $/;
  }

=head
Blast hit is at End of Query - 4 possibilites
                   valid     qs   qe    ss  se
start subject       yes      986  1000  1   15
start subject RC    no       986  1000  15  1
end   subject       no       986  1000  486 500
end   subject RC    yes      986  1000  500 486
=cut

  # End Q
  # 	extended seq will be Q->S
  #
  elsif ( $qe == $qlen && $ss == 1 ) { # Start S hits query    TESTED
    if ($subjseq->length() == $se) {   # occurs when Subject seq exists completely at end of Query seq
      #print STDERR "Subject $s exists completely at 3' end of $q\n";
      return undef;
    }
  
    $eseq = ">$q-$s overlap $qs $qe\n";
    $eseq .= $self->{queryseq} . substr($self->_getsubjseq($s)->seq, $se) . $/;
  }
  elsif ( $qe == $qlen && $ss == $slen) { # End S (RC) hits query      TESTED
    $eseq = ">$q-${s}rc overlap $qs $qe\n";
    $eseq .= $self->{queryseq} . substr($self->_getsubjseq($s)->revcom->seq, $ss - $se + 1) . $/;
  }

  return $eseq;
}


sub _getsubjseq {
  my ($self, $s) = @_;  
  #my $mm = `getfastaseqs.pl -f -s $self->{target} '$s'`;
  my $mm = `samtools faidx $self->{target} '$s'`;
  $mm =~ s/^>.+?\n//;
  $mm =~ s/\n//g;
  return Bio::Seq->new(-seq => $mm);
}






package main;
use strict;
use warnings;
use Getopt::Long;
#use ExtendSequence;


my %h;
my @keys = qw(query|q=s target|t=s blastdb|b=s num_threads|t=i minoverlap|m=i blastout|o=s pid|p=i help|h);
GetOptions (\%h, @keys);


die join(' ', @keys), $/ if $h{'help'};
if (!$h{'query'}) {
  die "Please define your 'query' sequence file\n";
}

my $es = ExtendSequence->new(\%h);
$es->run;
$es->report;

