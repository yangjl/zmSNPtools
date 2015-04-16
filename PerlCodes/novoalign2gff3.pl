#!/usr/bin/perl -w

# File: novoalign2gff3.pl
# Author: Eddy Yeh
#
# Description: This script takes as input the default alignment output file produced by Novoalign and
# generates a GFF3 file based on "unique" hits reported by novoalign.
#
# Change Log v1.2 (2010.04.28)
#  - Implemented parsing of variable mismatch numbers per read according to their
#    read length
#
# Change Log v1.1 (2010.03.12)
#  - Changed how mismatches are computed. Now, it actually counts every insertion base
#    as a mismatch (i.e. 8+GATCTC is no longer 1 mismatch, it's actually 6 mismatches)
#  - Fixed parsing when novoalign does not have quality values
#
# Change Log v1.0 (2009.11.07)
#  - Initial implementation of the script
#

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use POSIX qw(floor ceil);
use Term::ANSIColor;

use constant VERSION => "1.2 (2010.04.28)";
use constant true => 1;
use constant false => 0;
use constant SOURCE => "novoalign";
use constant FEATURE => "unique_hit";
use constant MAX_MISMATCH => 4;			# Default value of number of mismatches per read to allow

my ($file, @mismatches, $source, $feature, $help);
my $result = &GetOptions("help|h" => \$help,
                         "input|i=s" => \$file,
                         "mismatches|mismatch|mm:i{1,2}" => \@mismatches,
						 "source|s:s" => \$source,
						 "feature|f:s" => \$feature);

unless ($result && defined($file)) {
	print STDERR sprintf("\n");
	print STDERR colored("*** INCORRECT NUMBER OF ARGUMENTS ***", "bold red ON_BLACK") . "\n" if (!defined($help));
	print STDERR colored("USAGE:", "bold green ON_BLACK") . "\n";
	print STDERR colored(sprintf("   perl %s --input|-i <novoalign output file> [OPTIONS]", $0), "bold cyan ON_BLACK") . "\n\n";
	print STDERR colored("WHERE:", "bold green ON_BLACK") . "\n";
	print STDERR sprintf("   --input|-i <novoalign output file>    : Path to the alignment output file produced by novoalign\n");
	print STDERR sprintf("\n");
	print STDERR colored("OPTIONS:", "bold green ON_BLACK") . "\n";
	print STDERR sprintf("   --mismatches|-mm <number>             : Allow at most 'number' of mismatches uniformly across all reads\n");
	print STDERR sprintf("                                           regardless of read length. Insertions, deletions, and\n");
	print STDERR sprintf("                                           substitutions are counted as mismatches. This is the default\n");
	print STDERR sprintf("                                           parsing option [DEFAULT: %d]\n", MAX_MISMATCH);
	print STDERR sprintf("   --mismatches|-mm <number> <length>    : Allow at most 'number' of mismatches per 'length' of read. Insertions,\n");
	print STDERR sprintf("                                           deletions, and substitutions are counted as mismatches. This option\n");
	print STDERR sprintf("                                           allows parsing variable number of mismatches of reads with different\n");
	print STDERR sprintf("                                           lengths in the alignment output file. For instance, if '--mismatches 2 32'\n");
	print STDERR sprintf("                                           is specified, this means that we are allowing 2 mismatches per every 32 bp\n");
	print STDERR sprintf("                                           of the read. If the read is 75 bp in length, the total number of mismatches\n");
	print STDERR sprintf("                                           allowed is: CEILING((75 x 2) / 32) = 5 mismatches\n");
	print STDERR sprintf("   --source|-s <string>                  : Source string to include in GFF3 output (DEFAULT: %s)\n", SOURCE);
	print STDERR sprintf("   --feature|-f <string>                 : Feature string to include in GFF3 output (DEFAULT: %s)\n", FEATURE);
	print STDERR sprintf("   --help| -h                            : Display usage and command line arguments\n");
	print STDERR sprintf("\n");
	print STDERR colored("NOTE:", "bold green ON_BLACK") . " Only " .  colored("\"unique\"", "bold yellow ON_BLACK");
	print STDERR sprintf(" hits reported by novoalign are considered.\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("VERSION: %s\n", VERSION);
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

# Assigning default values
if (scalar(@mismatches) == 0) {
	push(@mismatches, MAX_MISMATCH);
} # End of if statemnet
$source = SOURCE if (!defined($source) || length($source) == 0);
$feature = FEATURE if (!defined($feature) || length($feature) == 0);

my $fh = new FileHandle();
open ($fh, $file) or die("Cannot open input file $file\n");

printf("##gff-version\t3\n");
printf("# %s\n", scalar(localtime(time)));
printf("#\n");
printf("# Parser: %s\n", $0);
printf("# Version: %s\n", VERSION);
printf("#\n");

if (scalar(@mismatches) == 1) {
	printf("# Uniform mismatches allowed per read: %s (regardless of read length)\n", $mismatches[0]);
} # End of if statement
else {
	printf("# Non-Uniform mismatches allowed per read: %s %s (allow at most %s mismatches per every %s bp)\n",
	       $mismatches[0], $mismatches[1], $mismatches[0], $mismatches[1]);
} # End of else statement

printf("#\n");
printf("# NOTE: Only unique hits (status \"U\") from default novoalign output are considered\n");
printf("\n");

while (<$fh>) {
	chomp;
	if (length($_) != 0 && $_ !~ m/^#/) {
		my ($read, $readtype, $readseq, $readqual, $status, $score, $qual, $ref, $refstart, $strand, $pseq, $poffset, $pstrand, $poly) = split(/\t/, $_);

		if ($status =~ m/^U$/i) {		# Only consider unique hits
			my $refend = $refstart + length($readseq) - 1;		# Need to adjust end position based on insertions

			# Fixing name if starts with ">"
			$read = $1 if ($read =~ m/^>(\S+)/i);
			$ref = $1 if ($ref =~ m/^>(\S+)/i);

			# Fix strand value
			$strand = $strand =~ m/^F$/i ? "+" : "-";

			my ($mis_count, $ins_count, $del_count) = (0, 0, 0);

			# Counting number of mismatches
			if (defined($poly) && length($poly) != 0) {
				foreach my $p (split(/\s/, $poly)) {
					if ($p =~ m/^(\d+)(.)>(.)$/i) {	# mismatch polymorphism
						$mis_count++;	# Increment by 1
					} # End of if statemnet
					elsif ($p =~ m/^(\d+)\+(\S+)$/i) {	# Insertion
						$ins_count += length($2);

						# only adjust refend if the insertion is not at the end
						# of the read sequence to avoid double counting
						$refend += length($2)if (length($readseq) - length($2) != $1);
					} # End of elsif statement
					elsif ($p =~ m/^(\d+)-(\S+)$/) {
						$del_count += length($2);
						# adjust refend based on number of deleted nucleotides
						$refend -= length($2);
					} # End of else statement
				} # End of for each statement
			} # End of if statement

			# validate number of mismatches
			my $mis_allowed = scalar(@mismatches) == 1 ? $mismatches[0] : ceil( (length($readseq) * $mismatches[0]) / $mismatches[1]);
			if ($mis_count + $ins_count + $del_count <= $mis_allowed) {	# OK to print
				printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tName=%s;Note=Polymorphisms: %d (%d mismatches, %s insertions, %s deletions)\n",
				       $ref, $source, $feature, $refstart, $refend, ".", $strand, ".",
					   $read, $mis_count+$ins_count+$del_count, $mis_count, $ins_count, $del_count);
			} # End of if statement
		} # End of if statement
	} # End of if statement
} # end of while loop
close ($fh);

