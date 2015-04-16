#!/usr/bin/perl -w

# Copyright (c) 2010-2011
# Cheng-Ting "Eddy" Yeh (eddyyeh@iastate.edu)
# Schnable Laboratory
# Iowa State University
# All Rights Reserved
#
# This script takes a fastq variant as input and converts the quality values to a different fastq variant
#
# Change Log v0.01 - 2011.03.28 (Eddy)
#  o Initial implementation of this script

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil floor);
use FileHandle;

use constant VERSION => "0.01 (2011.03.28)";
use constant true => 1;
use constant false => 0;
use constant DEFAULT_PROCESSORS_COUNT => 1;
use constant DEFAULT_BLOCK_SIZE => 700000;

# remove leading and trailing white space characters from string
sub trim {
    my $str = $_[0];
    
    if (defined($str)) {
        $str =~ s/^(\s+|\t+)//g;
        $str =~ s/(\s+|\t+)$//g;
    } # End of if statement

    return $str;
} # End of sub trim

# Rounds the decimal number to integer
sub round {
    my $number = shift;
    return int($number + 0.5);
}

# Converts fastq-sanger qualities to phred qualities
sub fastq_sanger_to_phred {
    my @qual = split(//, $_[0]);
    @qual = map { ord($_) - 33 } @qual;
    
    return @qual;
} # End of sub fastq_sanger_to_phred

# Converts phred qualities to fastq-sanger
sub phred_to_fastq_sanger {
    my @qual = @_; 
    @qual = map { chr($_ + 33) } @qual;

    return join("", @qual);
} # end of phred_to_fastq_sanger

# Converts fastq-solexa (prior 1.3) to phred
sub fastq_solexa_to_phred {
    my @qual = split(//, $_[0]);
    @qual = map{ int((10 * log(1 + 10 ** ((ord($_) - 64) / 10.0))) / log(10)) } @qual;

	return @qual;
} # End of fastq_solexa_to_phred

# Converts phred qualities to fastq-solexa (prior 1.3)
sub phred_to_fastq_solexa {
    my @qual = @_;
    @qual = map { chr(ceil(10 * (log(exp((log(10) * $_) / 10.0) - 1) / log(10)) + 64)) } @qual;

    return join("", @qual);
} # end of phred_to_fastq_solexa

# Converts fastq-illumina (1.3+) to phred qualities
sub fastq_illumina_to_phred {
    my @qual = split(//, $_[0]);
    @qual = map { ord($_) - 64 } @qual;

    return @qual;
} # End of sub fastq_illumina_to_phred

# Converts phred qualities to fastq-illumina (1.3+)
sub phred_to_fastq_illumina {
    my @qual = @_; 
    @qual = map { chr($_ + 64) } @qual;

    return join("", @qual);
} # End of sub phred_to_fastq_illumina

# Given a fastq file handler, retrieves the next available sequene and qualities
sub nextFastqSequence {
    my $fh = $_[0];
    my ($seq_desc, $seq, $qual_desc, $qual);
    my ($line, $stop, $prev_pos);

    # skip all the lines from input file handle until it reaches a line
    # whose first character is '@' denoting the sequence descriptor line
    $stop = false;
    while (!$stop && !eof($fh)) {
        $line = <$fh>;
        chomp($line);

        if ($line =~ m/^@/) {   # found sequence descritor
            $seq_desc = $line;
            $stop = true;
        } # End of if statement
    } # End of while loop

    # the following lines are sequences, read all of them until a line
    # whose first character is '+' denoting the quality descritor, adjust
    # file position accordingly
    $prev_pos = tell($fh);  # get current position
    $stop = false;
    while (!$stop && !eof($fh)) {
        $line = <$fh>;
        chomp($line);

        if (length($line) != 0) {
            if ($line =~ m/^\+/) {  # found quality descritor
                # adjust line accordinly
                seek($fh, $prev_pos, 0); 
                $stop = true;
            } # End of if statement
            elsif (!defined($seq) || length($seq) == 0) {
                $seq = $line;
            } # End of elsif statement
            else {
                $seq .= $line;
            } # end of else statement
        } # End of if statement

        $prev_pos = tell($fh);
    } # end of while loop

    # skip all the lines from input file handle until it reaches a line
    # whose first character is '+' denoting the quality descriptor line
    $stop = false;
    while (!$stop && !eof($fh)) {
        $line = <$fh>;
        chomp($line);

        if ($line =~ m/^\+/) {  # found quality descriptor
            $qual_desc = $line;
            $stop = true;
        } # end of if statement
    } # End of while loop

    # the following lines are qualities, read all of them until a line
    # whose first character is '@' is encountered or end of file is
    # encountered if and only if length of quality is equal or 
    # greater than length of sequence
    $prev_pos = tell($fh);
    $stop = false;
    while (!$stop && !eof($fh)) {
        $line = <$fh>;
        chomp($line);

        if (length($line) != 0) {
            if ($line =~ m/^@/ && defined($qual) && length($qual) >= length($seq)) {
                # adjust line accordinly
                seek($fh, $prev_pos, 0);
                $stop = true;
            } # End of if statement
            elsif (!defined($qual) || length($qual) == 0) {
                $qual = $line;
            } # End of else if statement
            else {
                $qual .= $line;
            } # End of else statement
        } # End of if statement

        $prev_pos = tell($fh);
    } # End of while loop

	my $sequence;
	# Got all the pieces we need for fastq sequence, check if they are correct
	if (defined($seq_desc) && defined($seq) && defined($qual_desc) && defined($qual)) {
		if (length($qual_desc) > 1 && substr($qual_desc, 1) ne substr($seq_desc, 1)) {
			print STDERR sprintf("\n");
			print STDERR sprintf("  ERROR: Invalid sequence and quality descriptor in fastq file\n");
			print STDERR sprintf("      o SEQUENCE DESCRIPTOR : %s\n", $seq_desc);
			print STDERR sprintf("      o QUALITY DESCRIPTOR  : %s\n", $qual_desc);
			print STDERR sprintf("\n");
			exit();
		} # End of if statemnet

		if (length($seq) != length($qual)) {
			print STDERR sprintf("\n");
			print STDERR sprintf("   ERROR: Length of sequence and quality does not match\n");
			print STDERR sprintf("      o DESCRIPTOR : %s\n", substr($seq_desc, 1));
			print STDERR sprintf("      o SEQUENCE   : %s (%s)\n", $seq, length($seq));
			print STDERR sprintf("      o QUALITY    : %s (%s)\n", $qual, length($qual));
			print STDERR sprintf("\n");
			exit();
		} # End of if statement
	} # End of if statemnet

	$seq_desc = substr($seq_desc, 1);
	$qual_desc = substr($qual_desc, 1);

    return ($seq_desc, $seq, $qual_desc, $qual);
} # end of sub nextFastqSequence

# prints the specified fastq pieces accordingly
sub formatFastq {
    my ($handle, $seq_name, $seq, $qual_name, $qual) = @_;  
    
	print $handle sprintf("@%s\n%s\n+%s\n%s\n", $seq_name, $seq, $qual_name, $qual);
} # End of sub formatFastq 

# Converts given fastq variant
sub convertFastq {
	my ($in_file, $in_qual_function, $out_file, $out_qual_function) = @_;

	# Open fastq file handler
	my $inFH = new FileHandle();
	open ($inFH, $in_file) or die("Cannot open fastq file '$in_file'\n");

	# Open fasta file handler
	my $outFH = new FileHandle();
	open ($outFH, sprintf(">%s", $out_file)) or die("Cannot create fasta file '$out_file'\n");

	my ($seq_desc, $seq, $qual_desc, $qual);
	my @qualities;
	while (!eof($inFH)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($inFH);
		@qualities = &$in_qual_function($qual);
		$qual = &$out_qual_function(@qualities);

		&formatFastq($outFH, $seq_desc, $seq, $qual_desc, $qual);
	} # End of while loop

	close ($inFH);
	close ($outFH);
} # end of sub convertFastq

my ($in_file, $it, $out_file, $ot, $processors, $block);

my $result = &GetOptions("in|i=s{1}" => \$in_file,
                         "it=s{1}" => \$it,
					     "out|o=s{1}" => \$out_file,
					     "ot=s{1}" => \$ot,
						 "processors|p:i{1}" => \$processors,
						 "block|b:i{1}" => \$block);

unless ($result && defined($in_file) && defined($out_file) &&
        $it =~ m/^(fastq-sanger|fastq-solexa|fastq-illumina)$/i &&
		$ot =~ m/^(fastq-sanger|fastq-solexa|fastq-illumina)$/i) {
	print STDERR sprintf("\n");
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("   perl %s -in <fastq file> -it <fastq variant> -out <fastq file> -ot <fastq-variant> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --in <fastq file>             : Input fastq file to be converted\n");
	print STDERR sprintf("   --it <type>                   : Specify the input fastq variant:\n");
	print STDERR sprintf("                                   'fastq-sanger'    : Sanger style fastq using PHRED and\n");
	print STDERR sprintf("                                                       ASCII offset of 33 (Range: 0 to 93)\n");
	print STDERR sprintf("                                   'fastq-solexa'    : Old solexa (prior to 1.3) quality values\n");
	print STDERR sprintf("                                                       with ASCII offset of 64 (Range: -5 to 63)\n");
	print STDERR sprintf("                                   'fastq-illumina'  : New Illumina/Solexa 1.3+ quality styles\n");
	print STDERR sprintf("                                                       with ASCII offset of 64 (Range: 0 to 63)\n");
	print STDERR sprintf("   --out <fastq file>            : Input fastq file to be converted\n");
	print STDERR sprintf("   --ot <type>                   : Specify the output fastq variant:\n");
	print STDERR sprintf("                                   'fastq-sanger'    : Sanger style fastq using PHRED and\n");
	print STDERR sprintf("                                                       ASCII offset of 33 (Range: 0 to 93)\n");
	print STDERR sprintf("                                   'fastq-solexa'    : Old solexa (prior to 1.3) quality values\n");
	print STDERR sprintf("                                                       with ASCII offset of 64 (Range: -5 to 63)\n");
	print STDERR sprintf("                                   'fastq-illumina'  : New Illumina/Solexa 1.3+ quality styles\n");
	print STDERR sprintf("                                                       with ASCII offset of 64 (Range: 0 to 63)\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("   --processors or -p <num>      : Specify the number of processors to use [DEFAULT: %s]\n", DEFAULT_PROCESSORS_COUNT);
	print STDERR sprintf("   --block or -b <size>          : Specify the number of reads to be read as blocks when splitting\n");
	print STDERR sprintf("                                   the input file for multiprocessor run [DEFAULT: %s]\n", DEFAULT_BLOCK_SIZE);
	print STDERR sprintf("\n");
	print STDERR sprintf("VERSION: %s\n", VERSION);
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

$it = lc($it);
$ot = lc($ot);
$processors = DEFAULT_PROCESSORS_COUNT if (!defined($processors) || $processors !~ m/^\d+$/);
$block = DEFAULT_BLOCK_SIZE if (!defined($block) || $block !~ m/^\d+$/);

if ($it eq $ot) {
	print STDERR sprintf("\n");
	print STDERR sprintf("Input fastq variant cannot be the same as the output fastq variant\n");
	print STDERR sprintf("\n");
	exit();
} # end of if statement

my ($inFH, $outFH);
my $id = $$;		# Process ID
my $input_qual_function;
my $output_qual_function;

# Determine the fastq-variant from input
if ($it eq "fastq-sanger") {
	$input_qual_function = \&fastq_sanger_to_phred;
} # end of if statement
elsif ($it eq "fastq-solexa") {
	$input_qual_function = \&fastq_solexa_to_phred;
} # end of else if statement
elsif ($it eq "fastq-illumina") {
	$input_qual_function = \&fastq_illumina_to_phred;
} # end of else if statement
else {
	die("Invalid fastq-variant specified for input\n");
} # end of else statement

# Determine the fastq-variant for output
if ($ot eq "fastq-sanger") {
	$output_qual_function = \&phred_to_fastq_sanger;
} # end of if statement
elsif ($ot eq "fastq-solexa") {
	$output_qual_function = \&phred_to_fastq_solexa;
} # end of else if statement
elsif ($ot eq "fastq-illumina") {
	$output_qual_function = \&phred_to_fastq_illumina;
} # end of else if statement
else {
	die("Invalid fastq-variant specified for output\n");
} # end of else statement

if ($processors == 1) {	# Uni-processor run
	&convertFastq($in_file, $input_qual_function, $out_file, $output_qual_function);
} # end of if statement
else {
	my %handlers;		# To keep track of different handlers
	my ($fh, $ofh);		# I/O Handllers
	my ($seq_desc, $seq, $qual_desc, $qual);

	$fh = new FileHandle();
	open ($fh, $in_file) or die("Cannot open input file '$in_file'\n");
	
	my $pnum = 0;
	my $count = 0;
	while (!eof($fh)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fh);

		if (!exists $handlers{$pnum}) {	# Create temp file
			$ofh = new FileHandle();
			open ($ofh, sprintf(">.%s.%s", $pnum, $id)) or die("Cannot create temporary file\n");
			$handlers{$pnum} = $ofh;
		} # end of if statement
		else {
			$ofh = $handlers{$pnum};
		} # end of else statement

		&formatFastq($ofh, $seq_desc, $seq, $qual_desc, $qual);
		$count++;

		if ($count % $block == 0) {
			$pnum = 0 if (++$pnum == $processors);	# Reset processor id if necessary
		} # End of if statement
	} # end of while loop
	close ($fh);

	my @pids = sort {$a <=> $b} keys %handlers;

	# Close all open handlers
	foreach my $p (@pids) {
		close ($handlers{$p});
	} # end of for each statement

	# Initialize multiprocessor run
	my @children;
	foreach my $p (@pids) {
		my $pid = fork();		# Spawn a new process

		if ($pid) {	# Parent process
			push(@children, $pid);
		} # end of if statemetn
		else {
			&convertFastq(sprintf(".%s.%s", $p, $id), $input_qual_function, sprintf(".%s.%s.out", $p, $id), $output_qual_function);
			exit(0);
		} # end of else statement
	} # end of for each statement

	# Waiting for all 'forked' processes to finish
	foreach my $child (@children) {
		waitpid($child, 0);
	} # end of for each statement
	
	# Open output handler
	my $outFH = new FileHandle();
	open ($outFH, sprintf(">%s", $out_file)) or die("Cannot create output file\n");

	# Merge fasta and qualities
	foreach my $p (@pids) {
		# Merge fasta
		$inFH = new FileHandle();
		open ($inFH, sprintf(".%s.%s.out", $p, $id)) or die("Cannot open temporary fasta file\n");
		while (<$inFH>) {
			chomp;
			print $outFH sprintf("%s\n", $_);
		} # end of while loop
		close ($inFH);

		unlink(sprintf(".%s.%s", $p, $id));
		unlink(sprintf(".%s.%s.out", $p, $id));
	} # end of for each statement

	close ($outFH);
} # end of else statement
