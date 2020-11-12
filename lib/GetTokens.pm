
=head1 NAME


=head1 SYNOPSIS


=head1 EXAMPLE


=head1 DESCRIPTION 


=head1 TODO


=head1 FEEDBACK



=head1 BUGS

  Contact Paul Cantalupo pcantalupo_at_gmail-dot-com
  
=head1 AUTHOR

  Paul Cantalupo
  
=head1 VERSIONS

  0.01

=cut

# Let the code begin
package GetTokens;
use strict;
use warnings;
use Exporter;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(GetTokens
                  );




sub GetTokens {
  my ($tokenfile) = @_;

  my @tokens;
  
  # get tokens from file
  if ($tokenfile && -s $tokenfile) {
     open (INFILE, "<", $tokenfile) or die ("Can't open $tokenfile: $!\n");
     while (<INFILE>) {
       next if ($_ =~ /^#/ || $_ =~ /^$/);
       chomp;
       push (@tokens, $_);
     }
     close (INFILE);
  }

  # treat remaining arguments on cmd line as tokens
  while (my $token = shift @ARGV) {
    push (@tokens, $token);
  }

  # grab tokens from STDIN (only works if STDIN comes from a pipe)
  if (-p STDIN) {
    while (<>) {
      chomp;
      push (@tokens, $_);
    }
  }
  
  return @tokens;
}



1;

