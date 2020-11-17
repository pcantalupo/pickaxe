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

