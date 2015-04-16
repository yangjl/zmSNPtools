#!/usr/bin/perl -w

# File: gmap2gff3.pl
# Author: Eddy Yeh
#
# Description: This program parses the alignment output produced by GMAP into GFF3 file format
#              considering identity, coverage, and tail lengths
#
# Change Log v2.6 (2012.07.30 -- Eddy)
#   - Bug fix caused by updated version of GSNAP causing program not be able to
#     deternine EST start and end positions
#
# Change Log v2.5 (2012.07.06 -- Eddy)
#   - Bug fix caused by updated version of GMAP causing the program not be able to
#     determine the query length (EST length) from the standard output
#
# Change Log v2.4 (2010.04.16 -- Eddy)
#   - Changed the --ident and --coverage options and combined --percentage options
#   - Implemented parsing of non-uniform mismatches per read based on their
#     read length
# 
# Change Log v2.3 (2010.04.15 -- Eddy)
#   - Modified how sequence mismatch is computed. Before, mismatches were computed
#     by taking the sum of "mismatches", "indels", and "unknowns" directly from
#     GMAP output. This calculation only considers mismatches occurring in the
#     aligned regions and ignores the sequence tails. For this reason, mismatches
#     in this version is computed as "sequence length" - "matches"
#
# Change Log v2.2 (2010.03.05 -- Eddy)
#   - Added --unique option. If this option is specified, only unique alignments
#     are parsed
#
# Change Log v2.1 (2010.02.24 -- Eddy)
#   - Modified mismatch parameter which is the sum of "mismatches", "indels",
#     and "unknowns" from GMAP output
#   - Modified how coverage is computed. Previously, it was computed from
#     the start/end of the mRNA alignment. In version v2.1, the coverage
#     reflects the total number of aligned nucleotides (matches) divided
#     by the total length of the mRNA
#
# Change Log v2.0 (2010.02.09 -- Eddy)
#   - This is a new implementation of gmap2gff.pl parser to convert from GMAP output to
#     GFF3 file format independent from v1.0+
#   - Added command line argument --best to only print the best possible alignment
#     for an Genomic/EST pair. The best possible alignment is determined using
#     highest coverage, highest identity, and lowest number of nucleotides as tails
#   - Added command line argument --intron to print intron regions in the genomic sequence
#   - Added command line argument --mismatches to parse output according to the mismatch
#     counts rather than percentage/coverage values. If mismatches option is specified,
#     overwrites identity/coverage/tail options

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use POSIX qw(floor ceil);
use Term::ANSIColor;

use constant VERSION => "2.6 (2012.07.30)";
use constant true => 1;
use constant false => 0;
use constant DEFAULT_MIN_IDENTITY => 0.95;		# Default minimum identity for each exon
use constant DEFAULT_MAX_TAIL => 0.05;			# Default maximum tail percentage
use constant DEFAULT_ABSOLUTE_TAIL => 60;		# Default maximum absolute tail in bp
use constant DEFAULT_MIN_COVERAGE => 0.90;		# Default minimum coverage for EST
use constant DEFAULT_SOURCE => "GMAP";			# Default source string to be used in GFF3 file

# returns the minimum value from an array of numbers
sub min {
	my @sorted = sort {$a <=> $b} @_;
	return shift(@sorted);
} # End of sub min

# returns the maximum value from an array of numbers
sub max {
	my @sorted = sort {$a <=> $b} @_;
	return pop(@sorted);
} # end of sub max

sub computeTail {
	my ($gen_start, $gen_end, $gen_length, $est_start, $est_end, $est_length, $strand) = @_;
	my $gen_tail_front = $gen_start - 1;
	my $gen_tail_end = $gen_length - $gen_end;
	my $est_tail_front = $est_start - 1;
	my $est_tail_end = $est_length - $est_end;
	my $tail = 0;

	if ($strand eq "+") {		# EST in the same direction with genomic
		$tail = &min($gen_tail_front, $est_tail_front) + &min($gen_tail_end, $est_tail_end);
	} # End of if statement
	else {
		$tail = &min($gen_tail_front, $est_tail_end) + &min($gen_tail_end, $est_tail_front);
	} # End of else statement

	return $tail;
} # End of sub computeTail

# Validates using mismatches criteria
sub validateMismatches {
	my ($paths, $mismatches) = @_;
	my $result;

	foreach my $p (keys %{ $paths }) {
		my $tail = &computeTail($paths->{$p}->{"GENOMIC_START"}, $paths->{$p}->{"GENOMIC_END"}, $paths->{$p}->{"GENOMIC_LENGTH"},
		                        $paths->{$p}->{"EST_START"}, $paths->{$p}->{"EST_END"}, $paths->{$p}->{"EST_LENGTH"},
								$paths->{$p}->{"STRAND"});
		# Update mismatches = substitutions + tail
		$paths->{$p}->{"MISMATCHES"} += $tail;
		
		# Check if total number of mismatches is less than allowed
		if ($paths->{$p}->{"MISMATCHES"} <= $mismatches) {	# OK
			my $identity = $paths->{$p}->{"EST_IDENTITY"};
			my $coverage = ($paths->{$p}->{"EST_END"} - $paths->{$p}->{"EST_START"} + 1) / $paths->{$p}->{"EST_LENGTH"};

			my $next_id = scalar(keys %{ $result->{$coverage}->{$identity}->{$tail} }) + 1;
			$result->{$coverage}->{$identity}->{$tail}->{$next_id} = sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s", 
				"genomic", $paths->{$p}->{"GENOMIC"},  $paths->{$p}->{"GENOMIC_START"}, $paths->{$p}->{"GENOMIC_END"},
				$identity, $coverage, $tail, $paths->{$p}->{"MATCHES"}, $paths->{$p}->{"MISMATCHES"}, 
				$paths->{$p}->{"STRAND"}, $paths->{$p}->{"EST_START"}, $paths->{$p}->{"EST_END"},
				$paths->{$p}->{"EST_LENGTH"});

			my @positions = sort {$a <=> $b} keys %{ $paths->{$p}->{"EXONS"} };
			for (my $i=0; $i < scalar(@positions); $i++) {
				my $pos = $positions[$i];
				$result->{$coverage}->{$identity}->{$tail}->{$next_id} .= sprintf(";%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",
					"exon", $paths->{$p}->{"GENOMIC"}, $paths->{$p}->{"EXONS"}->{$pos}->{"GENOMIC_START"},
					$paths->{$p}->{"EXONS"}->{$pos}->{"GENOMIC_END"}, $paths->{$p}->{"EXONS"}->{$pos}->{"IDENTITY"},
					"NA", "NA", "NA", "NA", $paths->{$p}->{"STRAND"}, $paths->{$p}->{"EXONS"}->{$pos}->{"EST_START"},
					$paths->{$p}->{"EXONS"}->{$pos}->{"EST_END"}, $paths->{$p}->{"EST_LENGTH"});

				# Create intron information
				if ($i != scalar(@positions) - 1) {
					my $intron_start = $paths->{$p}->{"EXONS"}->{$positions[$i]}->{"GENOMIC_END"} + 1;
					my $intron_end = $paths->{$p}->{"EXONS"}->{$positions[$i+1]}->{"GENOMIC_START"} - 1;
					$result->{$coverage}->{$identity}->{$tail}->{$next_id} .= sprintf(";i%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",
						"intron", $paths->{$p}->{"GENOMIC"}, $intron_start, $intron_end, "NA", "NA", "NA", "NA", "NA",
						$paths->{$p}->{"STRAND"}, "NA", "NA", "NA") if ($intron_end - $intron_start + 1 > 0);
				} # End of if statement
			} # end of for each loop
		} # End of if statement
	} # End of for each statement

	return $result;
} # End of sub validateMismatches

# Validates using percentage criterias
sub validatePercentages {
	my ($paths, $min_identity, $min_coverage, $max_tail, $abs_tail, $ignoretail) = @_;
	my $result;

	foreach my $p (keys %{ $paths }) {
		# check if all exons are at least >= $min_identity
		# and coverage at least >= $min_coverage
		my $identity = $paths->{$p}->{"EST_IDENTITY"};
		my $exon_identity = $paths->{$p}->{"MIN_EXON_IDENTITY"};
		my $coverage = ($paths->{$p}->{"EST_END"} - $paths->{$p}->{"EST_START"} + 1) / $paths->{$p}->{"EST_LENGTH"};

		if ($exon_identity >= $min_identity && $coverage >= $min_coverage) {	# Continue processing
			# Check if tails <= min(EST Length * $max_tail, $abs_tail)
			my $tail_allowed = &min(int($paths->{$p}->{"EST_LENGTH"} * $max_tail), $abs_tail);	# Maximum allowed tail
			my $tail = &computeTail($paths->{$p}->{"GENOMIC_START"}, $paths->{$p}->{"GENOMIC_END"}, $paths->{$p}->{"GENOMIC_LENGTH"},
			                        $paths->{$p}->{"EST_START"}, $paths->{$p}->{"EST_END"}, $paths->{$p}->{"EST_LENGTH"},
									$paths->{$p}->{"STRAND"});

			if ($ignoretail || $tail <= $tail_allowed) {	# Save this entry, pass all criterias
				my $next_id = scalar(keys %{ $result->{$coverage}->{$identity}->{$tail} }) + 1;
				$result->{$coverage}->{$identity}->{$tail}->{$next_id} = sprintf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s", 
					"genomic", $paths->{$p}->{"GENOMIC"},  $paths->{$p}->{"GENOMIC_START"}, $paths->{$p}->{"GENOMIC_END"},
					$identity, $coverage, $tail, $paths->{$p}->{"MATCHES"}, $paths->{$p}->{"MISMATCHES"}, 
					$paths->{$p}->{"STRAND"}, $paths->{$p}->{"EST_START"}, $paths->{$p}->{"EST_END"},
					$paths->{$p}->{"EST_LENGTH"});

				my @positions = sort {$a <=> $b} keys %{ $paths->{$p}->{"EXONS"} };
				for (my $i=0; $i < scalar(@positions); $i++) {
					my $pos = $positions[$i];
					$result->{$coverage}->{$identity}->{$tail}->{$next_id} .= sprintf(";%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",
						"exon", $paths->{$p}->{"GENOMIC"}, $paths->{$p}->{"EXONS"}->{$pos}->{"GENOMIC_START"},
						$paths->{$p}->{"EXONS"}->{$pos}->{"GENOMIC_END"}, $paths->{$p}->{"EXONS"}->{$pos}->{"IDENTITY"},
						"NA", "NA", "NA", "NA", $paths->{$p}->{"STRAND"}, $paths->{$p}->{"EXONS"}->{$pos}->{"EST_START"},
						$paths->{$p}->{"EXONS"}->{$pos}->{"EST_END"}, $paths->{$p}->{"EST_LENGTH"});

					# Create intron information
					if ($i != scalar(@positions) - 1) {
						my $intron_start = $paths->{$p}->{"EXONS"}->{$positions[$i]}->{"GENOMIC_END"} + 1;
						my $intron_end = $paths->{$p}->{"EXONS"}->{$positions[$i+1]}->{"GENOMIC_START"} - 1;
						$result->{$coverage}->{$identity}->{$tail}->{$next_id} .= sprintf(";%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",
							"intron", $paths->{$p}->{"GENOMIC"}, $intron_start, $intron_end, "NA", "NA", "NA", "NA", "NA",
							$paths->{$p}->{"STRAND"}, "NA", "NA", "NA") if ($intron_end - $intron_start + 1 > 0);
					} # End of if statement
				} # end of for each loop
			} # End of if statement
		} # End of if statement
	} # End of for each statement

	return $result;
} # End of sub validatePercentages

# Print the items that meet specified criterias
sub print {
	my ($name, $result, $source, $match_count, $best, $unique, $intron) = @_;

	if (defined($result)) {
		if ($best) {
			my @coverage = sort {$a <=> $b} keys %{ $result };
			my $c = shift(@coverage);
			my @identity = sort {$a <=> $b} keys %{ $result->{$c} };
			my $i = shift(@identity);
			my @tail = sort {$a <=> $b} keys %{ $result->{$c}->{$i} };
			my $t = shift(@tail);
			my @keys = sort {$a <=> $b} keys %{ $result->{$c}->{$i}->{$t} };
			
			if (scalar(@keys) == 1) {
				my $k = shift(@keys);
				foreach my $entry (split(/;/, $result->{$c}->{$i}->{$t}->{$k})) {
					my ($type, $genomic, $gen_start, $gen_end, $identity, $coverage, $tail, 
						$matches, $mismatches, $strand, $est_start, $est_end, $est_length) = split(/,/, $entry);
					if ($type eq "genomic") {
						printf("%s\t%s\t%s\t%s\t%s\t%2.3f\t%s\t%s\t",
							   $genomic, $source, "match", $gen_start, $gen_end, $identity, $strand, ".");
						printf("ID=%s_Match%d;Name=%s;Note=length: %d bp, coverage: %2.3f, tail: %d bp, Matches: %d bp, Mismatches: %d bp;Target=%s %d %d\n",
							   $source, ++$match_count, $name, $est_length, $coverage, $tail, $matches, $mismatches, $name,
							   $est_start, $est_end);
					} # End of if statement
					elsif ($type eq "exon") {
						printf("%s\t%s\t%s\t%s\t%s\t%2.3f\t%s\t%s\t",
							   $genomic, $source, "HSP", $gen_start, $gen_end, $identity, $strand, ".");
						printf("Parent=%s_Match%d;Target=%s %d %d\n",
							   $source, $match_count, $name, $est_start, $est_end);
					} # End of else statement
					elsif ($intron && $type eq "intron") {
						printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
							   $genomic, $source, "INTRON", $gen_start, $gen_end, ".", $strand, ".");
						printf("Parent=%s_Match%d\n", $source, $match_count);
					} # End of else if statement
				} # End of for each statement
			} # End of if statement
		} # End of if statement
		elsif ($unique) {
			my $paths = 0;
			foreach my $c (sort {$a <=> $b} keys %{ $result }) {
				foreach my $i (sort {$a <=> $b} keys %{ $result->{$c} }) {
					foreach my $t (sort {$a <=> $b} keys %{ $result->{$c}->{$i} }) {
						my @keys = sort {$a <=> $b} keys %{ $result->{$c}->{$i}->{$t} };
						foreach my $k (@keys) {
							foreach my $entry (split(/;/, $result->{$c}->{$i}->{$t}->{$k})) {
								my ($type, $genomic, $gen_start, $gen_end, $identity, $coverage, $tail, 
									$matches, $mismatches, $strand, $est_start, $est_end, $est_length) = split(/,/, $entry);
								$paths++ if ($type eq "genomic");
							} # End of for each statement
						} # End of for each statement
					} # End of for each statement
				} # End of for each statement
			} # End of for each statement
		
			if ($paths == 1) {
				my @coverage = sort {$a <=> $b} keys %{ $result };
				my $c = shift(@coverage);
				my @identity = sort {$a <=> $b} keys %{ $result->{$c} };
				my $i = shift(@identity);
				my @tail = sort {$a <=> $b} keys %{ $result->{$c}->{$i} };
				my $t = shift(@tail);
				my @keys = sort {$a <=> $b} keys %{ $result->{$c}->{$i}->{$t} };
				my $k = shift(@keys);
				foreach my $entry (split(/;/, $result->{$c}->{$i}->{$t}->{$k})) {
					my ($type, $genomic, $gen_start, $gen_end, $identity, $coverage, $tail, 
						$matches, $mismatches, $strand, $est_start, $est_end, $est_length) = split(/,/, $entry);
					if ($type eq "genomic") {
						printf("%s\t%s\t%s\t%s\t%s\t%2.3f\t%s\t%s\t",
							   $genomic, $source, "match", $gen_start, $gen_end, $identity, $strand, ".");
						printf("ID=%s_Match%d;Name=%s;Note=length: %d bp, coverage: %2.3f, tail: %d bp, Matches: %d bp, Mismatches: %d bp;Target=%s %d %d\n",
							   $source, ++$match_count, $name, $est_length, $coverage, $tail, $matches, $mismatches, $name,
							   $est_start, $est_end);
					} # End of if statement
					elsif ($type eq "exon") {
						printf("%s\t%s\t%s\t%s\t%s\t%2.3f\t%s\t%s\t",
							   $genomic, $source, "HSP", $gen_start, $gen_end, $identity, $strand, ".");
						printf("Parent=%s_Match%d;Target=%s %d %d\n",
							   $source, $match_count, $name, $est_start, $est_end);
					} # End of else statement
					elsif ($intron && $type eq "intron") {
						printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
							   $genomic, $source, "INTRON", $gen_start, $gen_end, ".", $strand, ".");
						printf("Parent=%s_Match%d\n", $source, $match_count);
					} # End of else if statement
				} # End of for each statement
			} # End of if statement
		} # End of else if statement
		else {
			foreach my $c (sort {$a <=> $b} keys %{ $result }) {
				foreach my $i (sort {$a <=> $b} keys %{ $result->{$c} }) {
					foreach my $t (sort {$a <=> $b} keys %{ $result->{$c}->{$i} }) {
						my @keys = sort {$a <=> $b} keys %{ $result->{$c}->{$i}->{$t} };
						foreach my $k (@keys) {
							foreach my $entry (split(/;/, $result->{$c}->{$i}->{$t}->{$k})) {
								my ($type, $genomic, $gen_start, $gen_end, $identity, $coverage, $tail, 
									$matches, $mismatches, $strand, $est_start, $est_end, $est_length) = split(/,/, $entry);
								if ($type eq "genomic") {
									printf("%s\t%s\t%s\t%s\t%s\t%2.3f\t%s\t%s\t",
										   $genomic, $source, "match", $gen_start, $gen_end, $identity, $strand, ".");
									printf("ID=%s_Match%d;Name=%s;Note=length: %d bp, coverage: %2.3f, tail: %d bp, Matches: %d bp, Mismatches: %d bp;Target=%s %d %d\n",
										   $source, ++$match_count, $name, $est_length, $coverage, $tail, $matches, $mismatches, $name,
										   $est_start, $est_end);
								} # End of if statement
								elsif ($type eq "exon") {
									printf("%s\t%s\t%s\t%s\t%s\t%2.3f\t%s\t%s\t",
										   $genomic, $source, "HSP", $gen_start, $gen_end, $identity, $strand, ".");
									printf("Parent=%s_Match%d;Target=%s %d %d\n",
										   $source, $match_count, $name, $est_start, $est_end);
								} # End of else statement
								elsif ($intron && $type eq "intron") {
									printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
										   $genomic, $source, "INTRON", $gen_start, $gen_end, ".", $strand, ".");
									printf("Parent=%s_Match%d\n", $source, $match_count);
								} # End of else if statement
							} # End of for each statement
						} # End of for each statement
					} # End of for each statement
				} # End of for each statement
			} # End of for each statement
		} # End of if statement
	} # End of if statement

	return $match_count;
} # End of sub print

# =============== MAIN PROGRAM STARTS HERE =============== #
my ($file, @mismatches, @percentages, $max_tail, $abs_tail, $ignoretail, $best, $unique, $intron, $source);
my $result = &GetOptions("input|i=s" => \$file,
                         "mismatches|mm:i{1,2}" => \@mismatches,
						 "pecentages|p:f{2}" => \@percentages,
						 "tail|t:f" => \$max_tail,
						 "absolutetail|at:i" => \$abs_tail,
						 "ignoretail|it!" => \$ignoretail,
						 "best|b!" => \$best,
						 "unique|u!" => \$unique,
						 "intron|int!" => \$intron,
						 "source|s:s" => \$source);

unless ($result && defined($file) && !(defined($best) && defined($unique))) {
	print STDERR sprintf("\n");
	print STDERR colored("*** INCORRECT NUMBER OF ARGUMENTS ***", "bold red ON_BLACK") . "\n";
	print STDERR colored("VERSION:", "bold green ON_BLACK") . " " . colored(sprintf("%s\n", VERSION), "bold cyan ON_BLACK");
	print STDERR colored("USAGE:", "bold green ON_BLACK") . " ";
	print STDERR colored(sprintf("perl %s --input|-i <gmap output file> [OPTIONS]", $0), "bold cyan ON_BLACK") . "\n";
	print STDERR colored("WHERE:", "bold green ON_BLACK") . "\n";
	print STDERR sprintf("   --input|-i <standard GMAP output file> : Complete path to the file produced by GMAP with -A\n");
	print STDERR sprintf("                                             argument option\n");
	print STDERR colored("OPTIONS:", "bold green ON_BLACK") . "\n";
	print STDERR sprintf("   --mismatches|-mm <n>                   : Maximum number of nucleotides allowed as mismatches for\n");
	print STDERR sprintf("                                            parsing. When this option is included, it overwrites\n");
	print STDERR sprintf("                                            --percentages, --tail, and --absolutetail options\n");
	print STDERR sprintf("                                            when evaluating alignments. Note, this option assumes at most\n");
	print STDERR sprintf("                                            \"n\" per read regardless of its length. For non-uniform\n");
	print STDERR sprintf("                                            number of mismatches, see --mismatches <n> <length>\n");
	print STDERR sprintf("   --mismatches|-mm <n> <length>          : Maximum number of nucleotides allowed as mismatches per \"length\"\n");
	print STDERR sprintf("                                            of the read. When this option is included, it overwrites\n");
	print STDERR sprintf("                                            --percentages, --tail, and --absolutetail options. For instance,\n");
	print STDERR sprintf("                                            if \"--mismatches 2 36\" is specified, this means to allow 2\n");
	print STDERR sprintf("                                            mismatches per every 36 bp of read length. If the read is 75 bp in\n");
	print STDERR sprintf("                                            length, the total number of mismatches allowed is:\n");
	print STDERR sprintf("                                               %s\n", colored("CEILING((75 x 2) / 36) = 5 mismatches", "bold yellow ON_BLACK"));
	print STDERR sprintf("   --percentages|-p <identity> <coverage> : Specify the minimum identity per segment and overall EST coverage\n");
	print STDERR sprintf("                                            allowed. This is the default parsing option. For instance, if\n");
	print STDERR sprintf("                                            \"--percentage 0.95 0.90\" is specified, this means to allow at least\n");
	print STDERR sprintf("                                            95%% identity per segment and overall EST coverage must be at least\n");
	print STDERR sprintf("                                            90%%. [DEFAULT: %2.2f %2.2f]\n", DEFAULT_MIN_IDENTITY, DEFAULT_MIN_COVERAGE);
	print STDERR sprintf("   --tail|-t <percentage>                 : Maximum percentage based on the length of the EST sequence to allow\n");
	print STDERR sprintf("                                            as tail. Fragments extending beyond the boundaries of the reference\n");
	print STDERR sprintf("                                            sequence is not considered as tail [DEFAULT: %2.2f]\n", DEFAULT_MAX_TAIL);
	print STDERR sprintf("   --absolutetail|-at <n>                 : Maximum number of bases to be used as absolute EST tails when\n");
	print STDERR sprintf("                                            parsing alignments [DEFAULT: %d]\n", DEFAULT_ABSOLUTE_TAIL);
	print STDERR sprintf("   --ignoretail|-it|--noignoretail|-noit  : If specified, tail criterias are ignored during parsing [DEFAULT: --noignoretail]\n");
	print STDERR sprintf("   --best|-b|--nobest|-nob                : If specified, only parse best possible alignment [DEFAULT: --nobest]\n");
	print STDERR sprintf("   --unique|-u|--nounique|-nou            : If specified, only parse unique aligments [DEFAULT: --nounique]\n");
	print STDERR sprintf("   --intron|-int|--nointron|-noint        : If specified, triggers the output of intron regions [DEFAULT: --nointron]\n");
	print STDERR sprintf("   --source|-s <source str>               : \"Source\" string to be used in the GFF3 file [DEFAULT: %s]\n", DEFAULT_SOURCE);
	print STDERR sprintf("\n");
	print STDERR colored("NOTES:", "bold green ON_BLACK") . " If " . colored("--ignoretails", "bold yellow ON_BLACK") . " ";
	print STDERR sprintf("is not specified, the total number of bases to allow as tails is the MINIMUM value\n");
	print STDERR sprintf("       between tail percentage (--tail) and absolute tail (--absolustetail) options.\n");
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

# Validating arguments
if (scalar(@mismatches) > 0 && scalar(@percentages) > 0) {
	print STDERR sprintf("\n");
	print STDERR sprintf("INVALID OPTIONS:\n");
	print STDERR sprintf("   --mismatches and --percentages options cannot be both specified\n");
	print STDERR sprintf("\n");
	exit();
} # End of if statement

if ($best && $unique) {
    print STDERR sprintf("\n");
    print STDERR sprintf("INVALID OPTIONS:\n");
    print STDERR sprintf("   --best and --unique options cannot be both specified\n");
    print STDERR sprintf("\n");
    exit();
} # End of if statement

# Assigning default values
$intron = false if (!defined($intron));     # default is to omit intron printing
$best = false if (!defined($best));         # default is to print all possible alignments unless --best is specified
$unique = false if (!defined($unique));     # default is to print all possible alignments unless --unique is specified
$ignoretail = false if (!defined($ignoretail)); # Default does not ignore tails when evaluating with percentage criterias
$max_tail = DEFAULT_MAX_TAIL if (!defined($max_tail) || $max_tail !~ m/^(\d+)$/i);	# Default maximum tail
$abs_tail = DEFAULT_ABSOLUTE_TAIL if (!defined($abs_tail) || $abs_tail !~ m/^(\d+)$/i);
$source = DEFAULT_SOURCE if (!defined($source) || length($source) == 0);
if (scalar(@percentages) == 0 && scalar(@mismatches) == 0) {	# Use default percentages criteria
	push(@percentages, DEFAULT_MIN_IDENTITY);
	push(@percentages, DEFAULT_MIN_COVERAGE);
} # End of if statement

my ($name, $pathnum, $paths);		# Shared variables
my $match_count = 0;

# Open file for processing
my $fh = new FileHandle();
open ($fh, $file) or die("Cannot open GMAP output file\n");

printf("##gff-version\t3\n");
printf("# %s\n", scalar(localtime(time)));
printf("#\n");
printf("# Parser: %s\n", $0);
printf("# Parser version: %s\n", VERSION);
printf("#\n");

if (scalar(@mismatches) == 1) {
	printf("# Uniform mismatches allowed per read: %s (regardless of read length)\n", $mismatches[0]);
} # End of if statemnet
elsif (scalar(@mismatches) == 2) {
	printf("# Non-Uniform mismatches allowed per read: %s %s (allow at most %s mismatches per every %s bp)\n",
	       $mismatches[0], $mismatches[1], $mismatches[0], $mismatches[1]);
} # End of else if statemnet
else {
	printf("# Minimum identity percentage per segment: %2.3f\n", $percentages[0]);
	printf("# Minimum coverage percentage: %2.3f\n", $percentages[1]);
	if (!$ignoretail) {
		printf("# Max. allowed tail: min(%2.2f%% of EST length, %d) bp\n", $max_tail * 100, $abs_tail);
	} # End of if statement
	else {
		printf("# Max. allowed tail: [Criteria ignored because of \"ignoretails\" was specified]\n"); 
	} # End of else statement
} # End of else statement
printf("# Introns parsing: %s\n", $intron ? "Yes" : "No");
printf("# Parse all possible aligments: %s\n", !$best && !$unique ? "Yes" : "No");
printf("# Best alignment parsing: %s\n", $best ? "Yes" : "No");
printf("# Unique aligment parsing: %s\n", $unique ? "Yes" : "No");
printf("\n");

while (<$fh>) {
	chomp;
	if ($_ =~ m/^>(\S+)/i) {
		if (defined($name) && defined($paths)) {
			my $result;

			if (scalar(@mismatches) == 1) {	# Uniform mismatches across all reads
				$result = &validateMismatches($paths, $mismatches[0]);	
			} # End of if statement
			elsif (scalar(@mismatches) == 2) {	# Non-Uniform mismatches depending of read length
				my $max_mismatches = ceil(($paths->{$pathnum}->{"EST_LENGTH"} * $mismatches[0]) / $mismatches[1]);
				$result = &validateMismatches($paths, $max_mismatches);
			} # End of else if statemnet
			else {
				$result = &validatePercentages($paths, $percentages[0], $percentages[1], $max_tail, $abs_tail, $ignoretail);
			} # End of else statement

			$match_count = &print($name, $result, $source, $match_count, $best, $unique, $intron);
		} # end of if statement

		$name = $1;
		undef($pathnum);
		undef($paths);		# Clear all previous paths
	} # End of if statement
	elsif ($_ =~ m/\s+Path\s+(\d+):\s+query\s+/i || $_ =~ m/\s+Alignment\s+for\s+path\s+(\d+):$/i) {	# obtain path number
		$pathnum = $1;
	} # End of else if statement
	elsif ($_ =~ m/\s+Coverage:\s+(-?\d+\.?\d*)\s+\(query\s+length:\s+(\d+)\s+bp\)$/i) {	# obtain EST length
		$paths->{$pathnum}->{"EST_LENGTH"} = $2;
	} # End of elsif statement
	elsif ($_ =~ m/\s+Percent\s+identity:\s+(\d+\.?\d*)\s+\((\d+)\s+matches,\s+(\d+)\s+mismatches,\s+(\d+)\s+indels,\s+(\d+)\s+unknowns\)/i) {	
		# obtain overall EST identity and match/mismatch bases
		$paths->{$pathnum}->{"EST_IDENTITY"} = $1 / 100;
		$paths->{$pathnum}->{"MATCHES"} = $2;
		$paths->{$pathnum}->{"MISMATCHES"} = $paths->{$pathnum}->{"EST_LENGTH"} - $paths->{$pathnum}->{"MATCHES"};	# mismatches
	} # End of elsif statement
	elsif ($_ =~ m/\s+Accessions:\s+\S+\s+\(out\s+of\s+(\d+)\s+bp\)/i) {	# obtain genomic length
		$paths->{$pathnum}->{"GENOMIC_LENGTH"} = $1;
	} # End of else statement
	elsif ($_ =~ m/\s+([+-])(\S+):-{0,}(\d+)-{1,}(\d+)\s+\((\d+)-(\d+)\)\s+(\d+)\%/i) {	# obtain exon information
		my $gen_start = &min($3, $4);
		my $gen_end = &max($3, $4);
		my $est_start = &min($5, $6);
		my $est_end = &max($5, $6);
		my $identity = $7 / 100;

		# Saving strand and genomic sequene information
		if (!exists $paths->{$pathnum}->{"GENOMIC"}) {
			$paths->{$pathnum}->{"STRAND"} = $1;						# Strand
			$paths->{$pathnum}->{"GENOMIC"} = $2;						# Genomic name position
			$paths->{$pathnum}->{"GENOMIC_START"} = $gen_start;			# Genomic start position
			$paths->{$pathnum}->{"GENOMIC_END"} = $gen_end;				# Genomic end position
			$paths->{$pathnum}->{"EST_START"} = $est_start;				# EST start position
			$paths->{$pathnum}->{"EST_END"} = $est_end;					# EST end position
			$paths->{$pathnum}->{"MIN_EXON_IDENTITY"} = $identity;		#  Min. Exon identity
		} # End of if statement
		else {		# Update start/end values
			$paths->{$pathnum}->{"GENOMIC_START"} = &min($paths->{$pathnum}->{"GENOMIC_START"}, $gen_start);
			$paths->{$pathnum}->{"GENOMIC_END"} = &max($paths->{$pathnum}->{"GENOMIC_END"}, $gen_end);
			$paths->{$pathnum}->{"EST_START"} = &min($paths->{$pathnum}->{"EST_START"}, $est_start);
			$paths->{$pathnum}->{"EST_END"} = &max($paths->{$pathnum}->{"EST_END"}, $est_end);
			$paths->{$pathnum}->{"MIN_EXON_IDENTITY"} = &min($paths->{$pathnum}->{"MIN_EXON_IDENTITY"}, $identity);
		} # End of else statment

		# Saving exon information
		$paths->{$pathnum}->{"EXONS"}->{$gen_start}->{"GENOMIC_START"} = $gen_start;
		$paths->{$pathnum}->{"EXONS"}->{$gen_start}->{"GENOMIC_END"} = $gen_end;
		$paths->{$pathnum}->{"EXONS"}->{$gen_start}->{"EST_START"} = $est_start;
		$paths->{$pathnum}->{"EXONS"}->{$gen_start}->{"EST_END"} = $est_end;
		$paths->{$pathnum}->{"EXONS"}->{$gen_start}->{"IDENTITY"} = $identity;
	} # End of if statement
} # end of while loop
close ($fh);

# very last one
if (defined($name) && defined($paths)) {
	my $result;

	if (scalar(@mismatches) == 1) {	# Uniform mismatches across all reads
		$result = &validateMismatches($paths, $mismatches[0]);	
	} # End of if statement
	elsif (scalar(@mismatches) == 2) {	# Non-Uniform mismatches depending of read length
		my $max_mismatches = ceil(($paths->{$pathnum}->{"EST_LENGTH"} * $mismatches[0]) / $mismatches[1]);
		$result = &validateMismatches($paths, $max_mismatches);
	} # End of else if statemnet
	else {
		$result = &validatePercentages($paths, $percentages[0], $percentages[1], $max_tail, $abs_tail, $ignoretail);
	} # End of else statement

	$match_count = &print($name, $result, $source, $match_count, $best, $unique, $intron);
} # end of if statement
