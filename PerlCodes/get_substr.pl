# !usr/bin/perl -w
# Jinliang Yang 11-24-2009

use strict;
use warnings;
use Tie::File;

my $filename = 'pseudo';
tie my @pseudo, 'Tie::File', $filename or die "Can't tie '$filename' $!";

# The output file is seq_ch4_169to178M
my $out= 'seq_ch4_169to178M';
open(OUTFILE, ">$out") or die "Could not open the outputfile $out\n";

for my $i (0..$#pseudo) {
if ( $pseudo[$i] =~ />chr4/ ) {
my ($seq) = substr($pseudo[$i + 1], 169000000,9000000);
print OUTFILE ">chr4_169M-178M\n$seq\n";
}
}
close(OUTFILE);
exit;
