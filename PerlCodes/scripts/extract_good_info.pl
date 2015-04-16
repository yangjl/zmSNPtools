#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Schnablelab::Tools;

if (scalar(@ARGV) != 1) {
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("perl $0 <lucy sequence file>\n");
	exit();
} # End of else statement

my ($name, $left, $right, $seq);
my $fh = new FileHandle();
open ($fh, $ARGV[0]) or die("Cannot open lucy sequences file");
printf("# Seq. Name: Sequence name identifier\n");
printf("# Length (raw): Total length of the sequence before trimming\n");
printf("# Left Clipping: The start position of the \"trimmed\" sequence\n");
printf("# Right Clipping: The end position of the \"trimmed\" sequence\n");
printf("# Trim. Seq. Length: Length of trimmed sequence\n");
printf("#\n");
printf("# Seq. Name\tLength (raw)\tLeft Clipping\tRight Clipping\tTrim. Seq. Length\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>/) {
			if ($_ =~ m/^>(\S+)\s\d+\s\d+\s\d+\s(\d+)\s(\d+)$/) {
				if (defined($name)) {
					printf("%s\t%s\t%s\t%s\t%s\n", $name, length($seq), $left, $right, $right - $left + 1);
				} # End of if statement

				$name = $1;
				$left = $2;
				$right = $3;
				$seq = "";
			} # End of if statement
			else {
				undef($name);
				undef($left);
				undef($right);
				$seq = "";
			
				print STDERR sprintf("*** Invalid format '%s'\n", $_);
			} # End of if statement
		} # End of if statement
		else {
			$seq .= trim($_);
		} # end of else statement
	} # End of if statement
} # End of while loop
close ($fh);

if (defined($name)) {
	printf("%s\t%s\t%s\t%s\t%s\n", $name, length($seq), $left, $right, $right - $left + 1);
} # End of if statement

