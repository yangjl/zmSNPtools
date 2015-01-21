#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

my ($breakpoints, $fasta);

my $result = &GetOptions("breakpoints|b=s{1}" => \$breakpoints,
                         "fasta|f=s{1}" => \$fasta);

unless ($result && defined($breakpoints) && defined($fasta)) {
	print STDERR sprintf("perl %s --breakpoints <breakpoints file> --fasta <fasta>\n", $0);
	exit();
} # end of unless statemetn

# Reading breakpoints
my %contigs;
my $fh = new FileHandle();
open ($fh, $breakpoints) or die("Cannot open breakpoints file\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0 && $_ !~ m/^#/) {
		my ($id, $length, $source, $ref, $ref_start, $ref_end, $ident, $cvg, $tail_f, $tail_e, $overhang, $aligned,
		    $align_start, $align_end, $break_left, $break_right, $cat, $total_taol, $tail_percentage) = split(/\t/, $_);
			$contigs{$id}->{"category"} = $cat;
			$contigs{$id}->{"align_start"} = $align_start;
			$contigs{$id}->{"align_end"} = $align_end;
	} # end of if statement
} # End of while loop
close ($fh);

$fh = new FileHandle();
open ($fh, $fasta) or die("Cannot open fasta\n");
my ($id, $seq);
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($id) && !exists $contigs{$id}) { # The whole sequence wasn't aligned
				&formatSequence(\*STDOUT, $id, $seq);
			} # end of if staement
			else {
				if (defined($id) && $contigs{$id}->{"category"} !~ m/^NO_TAIL$/i) {
					my $start = $contigs{$id}->{"align_start"};
					my $end = $contigs{$id}->{"align_end"};
					my $length = $end - $start + 1;
					my $mask = "X" x $length;
					substr($seq, $start - 1, $length, $mask);
					&formatSequence(\*STDOUT, $id, $seq);
				} # end of if statement
			} # End of if statement

			$id = $1;
			$seq = "";
		} # end of if statement
		else {
			$seq .= $_;	
		} # end of else statement
	} # end of if statement
} # end of while loop
close ($fh);

if (defined($id) && !exists $contigs{$id}) { # The whole sequence wasn't aligned
	&formatSequence(\*STDOUT, $id, $seq);
} # end of if staement
else {
	if (defined($id) && $contigs{$id}->{"category"} !~ m/^NO_TAIL$/i) {
		my $start = $contigs{$id}->{"align_start"};
		my $end = $contigs{$id}->{"align_end"};
		my $length = $end - $start + 1;
		my $mask = "X" x $length;
		substr($seq, $start - 1, $length, $mask);
		&formatSequence(\*STDOUT, $id, $seq);
	} # end of if statement
} # End of if statement


