#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Schnablelab::Tools;

sub truncateSequence {
	my ($seq, $left, $right) = @_;
	$seq = substr($seq, $left - 1, $right - $left + 1);
	return $seq;
} # End of sub truncateSequence 

if (scalar(@ARGV) != 1) {
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("perl $0 <lucy sequence file>\n");
	exit();
} # End of else statement

my ($name, $left, $right, $seq);
my $fh = new FileHandle();
open ($fh, $ARGV[0]) or die("Cannot open lucy sequences file");
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>/) {
			if ($_ =~ m/^>(\S+)\s\d+\s\d+\s\d+\s(\d+)\s(\d+)$/) {
				if (defined($name)) {
					formatSequence(\*STDOUT, $name, &truncateSequence($seq, $left, $right));
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
	formatSequence(\*STDOUT, $name, &truncateSequence($seq, $left, $right));
} # End of if statement

