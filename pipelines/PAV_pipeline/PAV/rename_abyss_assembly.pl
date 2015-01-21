#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

use constant DEFAULT_MINIMUM_LENGTH => 300;
use constant DEFAULT_N_COUNT => 3;

my ($input, $output, $min_length, $n);

my $result = &GetOptions("input|i=s{1}" => \$input,
                         "output|o=s{1}" => \$output,
						 "length|l:i{1}" => \$min_length,
						 "num|n:i{1}" => \$n);

unless ($result && defined($input) && defined($output)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("   perl %s --input <*-scaffolds.fa> --output <out.fas> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --input|-i <fasta file>       : Input assembled fasta file in which sequences\n");
	print STDERR sprintf("                                   will be renamed and splitted\n");
	print STDERR sprintf("   --output|-o <out file>        : Path to where renamed and splitted contigs will be saved\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("   --length|-l <bases>           : Minimum contig length to keep [DEFAULT: %s]\n", DEFAULT_MINIMUM_LENGTH);
	print STDERR sprintf("   --num|-n <n charas>           : Minimum number of N characters allowed before split is\n");
	print STDERR sprintf("                                   performed [DEFAULT: %s]\n", DEFAULT_N_COUNT);
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

$min_length = DEFAULT_MINIMUM_LENGTH if (!defined($min_length) || $min_length !~ m/^\d+$/);
$n = DEFAULT_N_COUNT if (!defined($n) || $n !~ m/^\d+$/);

my $prefix;
my @tmp = split(/\//, $input);
if ($tmp[$#tmp] =~ m/-(contigs|scaffolds)\.fa$/) {
	$prefix = $`;
} # end of if statement
else {
	$prefix = "assembly";
} # end of else statement

my $fh = new FileHandle();
open ($fh, $input) or die("Cannot open input file\n");

my $ofh = new FileHandle();
open ($ofh, sprintf(">%s", $output)) or die("Cannot create output file\n");

my ($name, $seq);
my $index = 0;
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($name)) {
				if ($n != 0) {
					my @tmp = split(/N{\Q$n\E,}/, $seq);

					for (my $i=0; $i < scalar(@tmp); $i++) {
						if (length($tmp[$i]) >= $min_length) {
							&formatSequence($ofh, sprintf("%s_contig%s", $prefix, ++$index), $tmp[$i]);
						} # end of if statement
					} # end of for loop
				} # end of if statement
				else {
					if (length($seq) >= $min_length) {
						&formatSequence($ofh, sprintf("%s_contig%s", $prefix, ++$index), $seq);
					} # end of if statement
				} # end of else statement
			} # end of if statemetn

			$name = $1;
			$seq = "";
		} # End of if statement
		else {
			$seq .= $_;
		} # end of else statement
	} # end of if statement
} # end of while loop

if (defined($name)) {
	if ($n != 0) {
		my @tmp = split(/N{\Q$n\E,}/, $seq);

		for (my $i=0; $i < scalar(@tmp); $i++) {
			if (length($tmp[$i]) >= $min_length) {
				&formatSequence($ofh, sprintf("%s_contig%s", $prefix, ++$index), $tmp[$i]);
			} # end of if statement
		} # end of for loop
	} # end of if statement
	else {
		if (length($seq) >= $min_length) {
			&formatSequence($ofh, sprintf("%s_contig%s", $prefix, ++$index), $seq);
		} # end of if statement
	} # end of else statement
} # end of if statemetn

close ($fh);
close ($ofh);
