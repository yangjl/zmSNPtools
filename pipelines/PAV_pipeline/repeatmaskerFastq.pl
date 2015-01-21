#!/usr/bin/perl -w

# Author: Eddy Yeh (eddyyeh@iastate.edu)
# Schnable Laboratory
# Iowa State University
#
# Repeat mask fastq sequences by aligning them against a repeat database using GSNAP
# Required Programs/Applications:
#   (1) gsnap
#   (2) gsnap2gff3.pl script

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Digest::MD5 qw(md5_hex);
use Time::Local;

use constant VERSION => "0.1 (2013.07.24)";
use constant true => 1;
use constant false => 0;
use constant DEFAULT_FASTQ_QUALITY_TYPE => "sanger";	# Default fastq quality type
use constant DEFAULT_REPEAT_MASK => 2;					# Default repeat masking character. 1 = 'X' ; 2 = 'N' ; 3 = lowercase (softmasking)
use constant DEFAULT_MISMATCHES_ALLOWED => 3;			# Default number of mismatches allowed over a specified length
use constant DEFAULT_MISMATCHES_LENGTH => 40;			# Default minimum sequence length to apply mismatch criteria
use constant DEFAULT_TAILS_ALLOWED => 5;				# Default maximum number of bases allowed as tails
use constant DEFAULT_TAILS_LENGTH => 75;				# Default minimum equence length to apply tail criteria
use constant DEFAULT_MINIMUM_REPEAT_LENGTH => 30;		# Default minimum repeat length
use constant DEFAULT_PROCESSORS_COUNT => 1;				# Maximum number of processors to be used for masking

# Format seconds notation into hours:mins:seconds
sub formatTime {
    my $totaltime = $_[0];
    my $str;

    $str = sprintf("%02d:", $totaltime / 3600);    # total hours
    $totaltime = $totaltime % 3600;
    $str .= sprintf("%02d:", $totaltime / 60);      # total minutes
    $totaltime = $totaltime % 60;
    $str .= sprintf("%02d", $totaltime);            # total sconds
    return $str;
} # End of sub formatTime

# Formats the parametric number and adds commas every thousand, millions, etc
sub formatNumber {
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;

    return $_;  
} # End of sub formatNumber

# checks the log file form GSNAP output and determines whether the alignment finished properly or not
sub checkGSNAP {
    my $log = $_[0];
    my $lastline = `tail -n 1 $log`;    # Get last line from log output
    my $ok = false;
    my $count = 0;

    if ($lastline =~ m/^Processed (\d+)/) {
        $count = $1; 
        $ok = true;
    } # end of if statement

    return ($ok, $count);
} # end of sub checkGSNAP

# Reads specified GFF3 file and extracts valid repeat regions
sub getRepeatRegions {
	my ($minLength, $gff3) = @_;
	my %repeats;

	my $fh = new FileHandle();
	open ($fh, $gff3) or die("Cannot open alignments gff3 file\n");
	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my @fields = split(/\t/, $_);
			if ($fields[2] =~ m/^HSP$/i) {
				my ($read, $start, $end);
				if ($fields[$#fields] =~ m/Target=(\S+) (\d+) (\d+)/) {
					$read = $1;
					$start = $2;
					$end = $3;

					if ($end - $start + 1 >= $minLength) {	# Keep
						if (exists $repeats{$read}) {
							$repeats{$read} .= sprintf(";%s-%s", $start, $end);
						} # end of if statemnt
						else {
							$repeats{$read} = sprintf("%s-%s", $start, $end);
						} # end of else statement
					} # End of if statement
				} # end of if statmenet
			} # end of if statement
		} # end of if statmenet
	} # end of while loop
	close ($fh);

	return \%repeats;
} # end of getRepeatRegions

# mask sequence with mask 1 (X)
sub maskSequence1 {
	my ($seq, $start, $end) = @_;
	my $length = $end - $start + 1;
	my $mask = "X" x $length;
	substr($seq, $start - 1, $length, $mask);
	return $seq;
} # end of sub maskSequence1

# mask sequence with mask 2 (N)
sub maskSequence2 {
	my ($seq, $start, $end) = @_;
	my $length = $end - $start + 1;
	my $mask = "N" x $length;
	substr($seq, $start - 1, $length, $mask);
	return $seq;
} # end of sub maskSequence2

# mask sequence with mask 3 (lowercase)
sub maskSequence3 {
	my ($seq, $start, $end) = @_;
	my $length = $end - $start + 1;
	my $mask = lc(substr($seq, $start - 1, $length));
	substr($seq, $start - 1, $length, $mask);
	return $seq;
} # End of maskSequence3

# count repeats with mask 1
sub countRepeats1 {
	my $seq = $_[0];
	my $count = $seq =~ tr/X//;
	return $count;
} # end of sub countRepeats1

# count repeats with mask 2
sub countRepeats2 {
	my $seq = $_[0];
	my $count = $seq =~ tr/N//;
	return $count;
} # end of sub countRepeats2

# count repeats with mask 3
sub countRepeats3 {
	my $seq = $_[0];
	my $count = $seq =~ tr/a-z//;
	return $count;
} # end of countRepeats3

# ==================== MAIN PROGRAM STARTS HERE ==================== #
my ($input, $output, $db, $dbpath, $mask, $proc);
my ($quality, $minLength, @mismatches, @tails);

my $result = &GetOptions("input|i=s{1}" => \$input,
                         "output|o=s{1}" => \$output,
						 "db|d=s{1}" => \$db,
						 "dbpath|x=s{1}" => \$dbpath,
						 "quality|q:s{1}" => \$quality,
						 "mask|r:i{1}" => \$mask,
						 "length|l:i{1}" => \$minLength,
						 "mismatches|m:i{1,2}" => \@mismatches,
						 "tails|t:i{1,2}" => \@tails,
						 "proc|p:i{1}" => \$proc);

unless ($result && defined($input) && defined($output) && defined($db) && defined($dbpath)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("   perl %s --input <fastq file> --output <out file> --db <gsnap db name> --dbpath <gsnap db path> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --input|-i <fastq file>         : Fastq file to be screened against repeats\n");
	print STDERR sprintf("   --output|o <out file>           : Output file to save masked sequences\n");
	print STDERR sprintf("   --db|-d <database name>         : Repeat database name assigned when running 'gmap_build'\n");
	print STDERR sprintf("   --dbpath|-x <database path>     : Repeat database path assigned when running 'gmap_build'\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("   --quality|-q <string>           : Fastq quality protocol. Allowed values are 'illumina' and 'sanger'\n");
	print STDERR sprintf("                                        'illumina' (ASCII 64-126) (equivalent to -J 64 -j -31)\n");
	print STDERR sprintf("                                        'sanger'   (ASCII 33-126) (equivalent to -J 33 -j 0)\n");
	print STDERR sprintf("                                     This option is required when aligning fastq files using gsnap [ DEFAULT: %s ]\n",
	                                                         DEFAULT_FASTQ_QUALITY_TYPE);
	print STDERR sprintf("   --mask|-r <number>              : Select masking characters to be used to for repetitive sequences [ DEFAULT: %s ]\n",
	                                                         DEFAULT_REPEAT_MASK);
	print STDERR sprintf("                                         1 = Mask repeat regions with 'X'\n");
	print STDERR sprintf("                                         2 = Mask repeat regions with 'N'\n");
	print STDERR sprintf("                                         3 = Mask repeat regions by converting them to lowercase\n");
	print STDERR sprintf("   --length|-l <min length>        : Minimum repeat length [ DEFAULT: %s ]\n", DEFAULT_MINIMUM_REPEAT_LENGTH);
    print STDERR sprintf("   --mismatches|-m <num> <length>  : Consider repeat sequences by allowing 'num' mismatches per 'length'\n");
	print STDERR sprintf("                                     of the read. This option regulates how the similarity of repeat sequences\n");
	print STDERR sprintf("                                     with the sequencing read. For instance: if \"--mismatches 3 40\" is specified,\n");
	print STDERR sprintf("                                     this means that we are allowing 3 mismatches for every 40 bp of read length.\n");
    print STDERR sprintf("                                     If read length is 100 bp, the total number of mismatches allowed is:\n");
	print STDERR sprintf("                                     CEILING((100 x 3) / 40) = 8 mismatches. When second argument 'length' is\n");
	print STDERR sprintf("                                     be omitted, the script will allow at most 'num' mismatches between the repeat\n");
	print STDERR sprintf("                                     sequence and read regardless of length. [ DEFAULT: %s %s ]\n", 
	                                                           DEFAULT_MISMATCHES_ALLOWED, DEFAULT_MISMATCHES_LENGTH);
    print STDERR sprintf("   --tail|-t <n> <length>          : Maximum number of nucleotides allowed as tails (front + end) when aligning repeats\n");
	print STDERR sprintf("                                     to the per read relative to read length. For instance if \"--tail 5 75\" is\n");
    print STDERR sprintf("                                     specified, this means that we are allowing 5 bases as tails per every 75 bp\n");
    print STDERR sprintf("                                     of the read. If read length is 100 bp, the total number of tails allowed\n");
    print STDERR sprintf("                                     is: CEILING((100 * 5) / 75) = 7 bases. When second argument 'length' is\n");
	print STDERR sprintf("                                     omitted, then script will allow at most 'num' bases as tails regardless\n");
	print STDERR sprintf("                                     of read length. [ DEFAULT: %s %s ]\n", DEFAULT_TAILS_ALLOWED, DEFAULT_TAILS_LENGTH);
	print STDERR sprintf("    --proc|p <int>                 : Total number of processors to use when masking [ DEFAULT: %s ]\n", DEFAULT_PROCESSORS_COUNT);
	print STDERR sprintf("\n");
	print STDERR sprintf("VERSION: %s\n", VERSION);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

# Assigning default values
$quality = DEFAULT_FASTQ_QUALITY_TYPE if (!defined($quality));
$mask = DEFAULT_REPEAT_MASK if (!defined($mask) || $mask !~ m/^\d+$/);
$minLength = DEFAULT_MINIMUM_REPEAT_LENGTH if (!defined($minLength) || $minLength !~ m/^\d+$/);
$proc = DEFAULT_PROCESSORS_COUNT if (!defined($proc) || $proc !~ m/^\d+$/);

if (scalar(@mismatches) == 0) {
	push(@mismatches, DEFAULT_MISMATCHES_ALLOWED, DEFAULT_MISMATCHES_LENGTH);
} # end of if statement

if (scalar(@tails) == 0) {
	push(@tails, DEFAULT_TAILS_ALLOWED, DEFAULT_TAILS_LENGTH);
} # end of if statement

# Validate quality type
$quality = lc($quality);	# Convert to lowercase
if ($quality !~ m/^(sanger|illumina)$/i) {
	print STDERR sprintf("\n");
	print STDERR sprintf("ERROR: Fastq quality protocol must be a 'illumina' or 'sanger'\n");
	print STDERR sprintf("         'illumina' (ASCII 64-126) (equivalent to -J 64 -j -31)\n");
	print STDERR sprintf("         'sanger'   (ASCII 33-126) (equivalent to -J 33 -j 0)\n");
	print STDERR sprintf("\n");
	exit();
} # end of if statement

# Validate mask
if ($mask < 0 || $mask > 3) {
	print STDERR sprintf("\n");
	print STDERR sprintf("ERROR: Repeat mask (--mask) must be an integer value between 1 and 3 inclusive\n");
	print STDERR sprintf("         1 = Mask repeat regions with 'X'\n");
	print STDERR sprintf("         2 = Mask repeat regions with 'N'\n");
	print STDERR sprintf("         3 = Mask repeat regions by converting them to lowercase (softmasking)\n");
	print STDERR sprintf("\n");
	exit();
} # end of if statement

my $runid = md5_hex($$);
my $gsnap = `which gsnap`;					chomp($gsnap);
my $gsnap2gff3 = `which gsnap2gff3.pl`;		chomp($gsnap2gff3);

if (length($gsnap) == 0 || length($gsnap2gff3) == 0) {
	print STDERR sprintf("\n");
	print STDERR sprintf("ERROR: Program cannot continue because it could nto detect the locations of\n");
	print STDERR sprintf("       'gsnap' and/or 'gsnap2gff3.pl' in your system. Pleaes ensure you incorporate\n");
	print STDERR sprintf("       the executable locations in your environment variable PATH\n");
	print STDERR sprintf("\n");
	print STDERR sprintf(" o gsnap can be downloaded from: http://research-pub.gene.com/gmap/\n");
	print STDERR sprintf(" o gsnap2gff3.pl is the custom GSNAP parser to filter out confident alignments\n");
	print STDERR sprintf("\n");
	exit();
} # end of if statement

# Running
my $logFH = new FileHandle();
open ($logFH, sprintf(">%s.log", $output)) or die("Cannot create log file\n");

my $start_time = timelocal(localtime(time));

my $command;
print $logFH sprintf("# %s\n", scalar(localtime(time)));
print $logFH sprintf("\n");
print STDERR sprintf("\n");
print $logFH sprintf(" o GSNAP Path: %s\n", $gsnap);
print STDERR sprintf(" o GSNAP Path: %s\n", $gsnap);
print $logFH sprintf(" o gsnap2gff3.pl Path: %s\n", $gsnap2gff3);
print STDERR sprintf(" o gsnap2gff3.pl Path: %s\n", $gsnap2gff3);
print $logFH sprintf(" o Database Name: %s\n", $db);
print STDERR sprintf(" o Database Name: %s\n", $db);
print $logFH sprintf(" o Databae Path: %s\n", $dbpath);
print STDERR sprintf(" o Database Path: %s\n", $dbpath);
print $logFH sprintf(" o Mismatches Allowed: %s\n", join(" ", @mismatches));
print STDERR sprintf(" o Mismatches Allowed: %s\n", join(" ", @mismatches));
print $logFH sprintf(" o Tails Allowed: %s\n", join(" ", @tails));
print STDERR sprintf(" o Tails Allowed: %s\n", join(" ", @tails));
print $logFH sprintf(" o Input Fastq: %s\n", $input);
print STDOUT sprintf(" o Input Fastq: %s\n", $input);
print $logFH sprintf(" o Output Fastq: %s\n", $output);
print STDOUT sprintf(" o Output Fastq: %s\n", $output);

if ($mask == 1) {
	print $logFH sprintf(" o Repeat Mask: %s (X)\n", $mask);
	print STDERR sprintf(" o Repeat Mask: %s (X)\n", $mask);
} # end of if statement
elsif ($mask == 2) {
	print $logFH sprintf(" o Repeat Mask: %s (N)\n", $mask);
	print STDERR sprintf(" o Repeat Mask: %s (N)\n", $mask);
} # end of elsif statement
else {
	print $logFH sprintf(" o Repeat Mask: %s (Lowercase)\n", $mask);
	print STDERR sprintf(" o Repeat Mask: %s (Lowercase)\n", $mask);
} # End of else statement

print $logFH sprintf(" o Mininum Repeat Length: %s\n", $minLength);
print STDERR sprintf(" o Mininum Repeat Length: %s\n", $minLength);
print $logFH sprintf(" o Total Processors: %s\n", $proc);
print STDERR sprintf(" o Total Processors: %s\n", $proc);

print $logFH sprintf("     + Running alignments ... ");
print STDERR sprintf("     + Running alignments ... ");

my $alignment_start_time = timelocal(localtime(time));
$command = sprintf("%s -d %s -D %s -B 5 -m 15 -i 2 -N 0 -t %s -n 5 --quality-protocol %s --nofails %s > %s.gsnap 2> %s.gsnap.log",
                   $gsnap, $db, $dbpath, $proc, $quality, $input, $runid, $runid);
system($command);
my $alignment_end_time = timelocal(localtime(time));

my ($ok, $count) = &checkGSNAP(sprintf("%s.gsnap.log", $runid));
if ($ok) {	# Alignment was successful
	print $logFH sprintf("DONE [ %s ]\n", &formatTime($alignment_end_time - $alignment_start_time));
	print STDERR sprintf("DONE [ %s ]\n", &formatTime($alignment_end_time - $alignment_start_time));
	
	print $logFH sprintf("     + Parsing alignments ... ");
	print STDERR sprintf("     + Parsing alignments ... ");
	my $parsing_start_time = timelocal(localtime(time));
	$command = sprintf("%s --format gff3 -i %s.gsnap -o %s.gff3 -m %s -t %s -p %s 2> %s.gff3.log",
	                   $gsnap2gff3, $runid, $runid, join(" ", @mismatches), join(" ", @tails), $proc, $runid);
	system($command);
	my $parsing_end_time = timelocal(localtime(time));
	print $logFH sprintf("DONE [ %s ]\n", &formatTime($parsing_end_time - $parsing_start_time));
	print STDERR sprintf("DONE [ %s ]\n", &formatTime($parsing_end_time - $parsing_start_time));

	print $logFH sprintf("     + Reading repeat regions ... ");
	print STDERR sprintf("     + Reading repeat regions ... ");
	my $reading_repeats_start_time = timelocal(localtime(time));
	my $repeats = &getRepeatRegions($minLength, sprintf("%s.gff3", $runid));
	my $reading_repeats_end_time = timelocal(localtime(time));
	print $logFH sprintf("DONE [ %s ]\n", &formatTime($reading_repeats_end_time - $reading_repeats_start_time));
	print STDERR sprintf("DONE [ %s ]\n", &formatTime($reading_repeats_end_time - $reading_repeats_start_time));

	# Determine masking function
	my $maskSequence;
	my $countRepeats;
	if ($mask == 1) {
		$maskSequence = \&maskSequence1;
		$countRepeats = \&countRepeats1;
	} # end of if statmeent
	elsif ($mask == 2) {
		$maskSequence = \&maskSequence2;
		$countRepeats = \&countRepeats2;
	} # end of else if statement
	else {
		$maskSequence = \&maskSequence3;
		$countRepeats = \&countRepeats3;
	} # end of else statement

	print $logFH sprintf("     + Masking sequences ... ");
	print STDERR sprintf("     + Masking sequences ... ");
	my $masking_start_time = timelocal(localtime(time));
	my $fh = new FileHandle();
	open ($fh, $input) or die("Cannot open file input file\n");

	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s", $output)) or die("Cannot create output file\n");

	my $totalSequences = 0;
	my $totalBases = 0;
	my $totalRepeatBases = 0;
	
	while (!eof($fh)) {
		my $seqdesc = <$fh>;		chomp($seqdesc);
		my $seq = <$fh>;			chomp($seq);
		my $qualdesc = <$fh>;		chomp($qualdesc);
		my $qual = <$fh>;			chomp($qual);
		$totalSequences++;

		# Obtain id
		my $id = substr($seqdesc, 1);
		if ($id =~ m/(\S+)\s/) {
			$id = $1;
		} # End of if statement
		
		if (exists $repeats->{$id}) {
			foreach my $t (split(/;/, $repeats->{$id})) {
				my ($start, $end) = split(/-/, $t);
				$seq = &$maskSequence($seq, $start, $end);	
			} # end of for each statement
		} # end of if statement
		
		my $sequenceLength = length($seq);
		my $repeatBases = &$countRepeats($seq);

		$totalBases += $sequenceLength;
		$totalRepeatBases += $repeatBases;

		# Print unmasked/masked sequence
		print $ofh sprintf("%s %s/%s=%2.1f%% Repetitive\n%s\n%s\n%s\n",
		                   $seqdesc, $repeatBases, $sequenceLength, ($repeatBases / $sequenceLength) * 100,
						   $seq, $qualdesc, $qual);
	} # end of while loop
	close ($fh);
	my $masking_end_time = timelocal(localtime(time));
	print $logFH sprintf("DONE [ %s ]\n", &formatTime($masking_end_time - $masking_start_time));
	print STDERR sprintf("DONE [ %s ]\n", &formatTime($masking_end_time - $masking_start_time));

	# Reporting
	print $logFH sprintf("     + Summary\n");
	print STDERR sprintf("     + Summary\n");
	print $logFH sprintf("         - Total Sequences: %s (Average Read Length: %s bp)\n", &formatNumber($totalSequences),
	                     $totalSequences != 0 ? int($totalBases / $totalSequences) : 0);
	print STDERR sprintf("         - Total Sequences: %s (Average Read Length: %s bp)\n", &formatNumber($totalSequences),
	                     $totalSequences != 0 ? int($totalBases / $totalSequences) : 0);
	print STDERR sprintf("         - Total Sequences: %s\n", &formatNumber($totalSequences));
	print $logFH sprintf("         - Total Base Pairs: %s\n", &formatNumber($totalBases));
	print STDERR sprintf("         - Total Base Pairs: %s\n", &formatNumber($totalBases));
	print $logFH sprintf("         - Total Repeat Base Pairs: %s (%2.1f%%)\n", &formatNumber($totalRepeatBases),
	                     $totalBases != 0 ? ($totalRepeatBases / $totalBases) * 100 : 0);
	print STDERR sprintf("         - Total Repeat Base Pairs: %s (%2.1f%%)\n", &formatNumber($totalRepeatBases),
	                     $totalBases != 0 ? ($totalRepeatBases / $totalBases) * 100 : 0);
	my $end_time = timelocal(localtime(time));
	print $logFH sprintf("     + Run-Time: %s\n", &formatTime($end_time - $start_time));
	print STDERR sprintf("     + Run-Time: %s\n", &formatTime($end_time - $start_time));

	# Clean up
	unlink(sprintf("%s.gsnap", $runid));
	unlink(sprintf("%s.gsnap.log", $runid));
	unlink(sprintf("%s.gff3", $runid));
	unlink(sprintf("%s.gff3.log", $runid));
} # end of if statement
else {
	print $logFH sprintf("ERROR\n");
	print $logFH sprintf("\n");
	print $logFH sprintf("There was an error performing gsnap alignments. Please ensure that the database name\n");
	print $logFH sprintf("and database path are in the proper locations specified\n");
	print $logFH sprintf("\n");
	print STDERR sprintf("ERROR\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("There was an error performing gsnap alignments. Please ensure that the database name\n");
	print STDERR sprintf("and database path are in the proper locations specified\n");
	print STDERR sprintf("\n");
	
	unlink(sprintf("%s.gsnap", $runid));
	unlink(sprintf("%s.gsnap.log", $runid));
} # end of else statmeent

print STDERR sprintf("\n");
close ($logFH);
