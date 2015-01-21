#!/usr/bin/perl -w
#
# Usage: perl fastaclean.pl -f <fasta files> [-l min_length -bp bp_per_line]
# ====================================================================================================

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

use constant BP_PER_LINE => 70;
use constant MIN_SEQ_LENGTH => 1;

# Prints the usage error message of this program
sub printUsage() {
	printf("\n");
	printf("*** INCORRECT NUMBER OF ARGUMENTS OR FILE DOES NTO EXISTS OR NOT READABLE ***\n");
	printf("USAGE: perl %s -f <fasta files> [OPTIONS]\n", $0);
	printf("OPTIONS:   --bp n     : Number of bases per line (DEFAULT: %d)\n", BP_PER_LINE);
	printf("           --length n : Minimum length of sequence in base pairs (DEFAULT: %d)\n", MIN_SEQ_LENGTH);
	printf("\n");
} # End of sub printUsage()

# ------------------------- MAIN PROGRAM STARTS HERE ------------------------- #
my (@files, $bp, $minlength);
my $result = &GetOptions("file|f=s{1,}" => \@files,
		                 "b|bp:i{1}" => \$bp,
		                 "l|length:i{1}" => \$minlength);

if ($result && scalar(@files) > 0) {
	$bp = BP_PER_LINE if (!defined($bp));
	$minlength = MIN_SEQ_LENGTH if (!defined($minlength));

	foreach my $f (@files) {
		my $name;            # The name of the current sequence
		my $seq;             # The accumulated bases of $seqName

		# open input file and start reading
		my $fh = new FileHandle;
		open ($fh, $f) or die("Error opening input file '$f'\n");

		while (<$fh>) {
			chomp;
			if (length($_) != 0) {
				if ($_ =~ m/^>(\S+)/) {
					if (defined($name) && length($seq) >= $minlength) {
						&formatSequence(\*STDOUT, $name, $seq, $bp);
					} # end of if statement

					$name = $1;
					$seq = "";
				} # end of if statement
				else {
					$seq .= uc($_);
				} # end of else statement
			} # end of if statemnet
		} # end of while loop
		
		if (defined($name) && length($seq) >= $minlength) {
			&formatSequence(\*STDOUT, $name, $seq, $bp);
		} # end of if statement

		close $fh;                # Closing input file
	} # end of for each statement
} # End of if statement
else {
	&printUsage();
} # End of else statement

# --------------------------- END OF MAIN PROGRAM --------------------------- #
