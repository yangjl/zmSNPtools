#!/usr/bin/perl -w

# Copyright (c) 2010-2011
# Cheng-Ting "Eddy" Yeh (eddyyeh@iastate.edu)
# Schnable Laboratory
# Iowa State University
# All Rights Reserved
#
# This scripts takes as command line argument a fastq file and its corresponding variant
# and converts the contents to standard fasta and quality files
#
# Change Log v0.03 - 2011.03.02 (Eddy)
#  o Added program parameter --processors|-p to enable conversion of file formats using multiple processors
#  o Program speed optimizations
#
# Change Log v0.02 - 2010.07.08 (Eddy)
#  o Changed how fastq-variant qualities are converted to phred scales
#
# Change Log v0.01 - 2010.06.11 (Eddy)
#  o Initial implementation of this script

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil floor);
use FileHandle;

use constant VERSION => "0.03 (2011.03.02)";
use constant true => 1;
use constant false => 0;
use constant QUAL_PER_LINE => 25;
use constant DEFAULT_PROCESSORS_COUNT => 1;
use constant DEFAULT_BLOCK_SIZE => 500000;

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
    my ($seq_name, $seq, $qual_name, $qual) = @_;  
    
    return sprintf("@%s\n%s\n+%s\n%s\n", $seq_name, $seq, $qual_name, $qual);
} # End of sub formatFastq 

# format the sequence to the specified file handler
sub formatSequence {
    my ($handle, $name, $seq) = @_; 

    # Removing unwanted leading/trailing white space characters
    $name = &trim($name);
    $seq = &trim($seq);
    
    if ($name =~ m/^>/) {
        print $handle sprintf("%s\n", $name);
    } # end of if statement
    else {
        print $handle sprintf(">%s\n", $name);
    } # end of else statement

	print $handle sprintf("%s\n", $seq);
} # End of sub formatSequence

# Format the quality values to the specified file handler
sub formatQuality {
    my ($handle, $name, @qual) = @_; 
    
    if ($name =~ m/^>/) {
        print $handle sprintf("%s\n", $name);
    } # End of if statemnet
    else {
        print $handle sprintf(">%s\n", $name);
    } # end of else statement

    while (scalar(@qual) > 0) {
        my @sub = splice(@qual, 0, QUAL_PER_LINE);
        @sub = map { sprintf("%2s", int($_)) } @sub;
        print $handle sprintf("%s\n", join(" ", @sub));
    } # End of while loop
} # end of sub formatQuality

# Converts the input file to fasta/quality format
sub convertFastqToFasta {
	my ($fastq_file, $fasta_file, $qual_file, $input_qual_function) = @_;

	# Open fastq file handler
	my $fastqFH = new FileHandle();
	open ($fastqFH, $fastq_file) or die("Cannot open fastq file '$fastq_file'\n");

	# Open fasta file handler
	my $fastaFH = new FileHandle();
	open ($fastaFH, sprintf(">%s", $fasta_file)) or die("Cannot create fasta file '$fasta_file'\n");

	# Open quality file handler
	my $qualFH = new FileHandle();
	open ($qualFH, sprintf(">%s", $qual_file)) or die("Cannot create quality file '$qual_file'\n");

	my ($seq_desc, $seq, $qual_desc, $qual);
	my @qualities;
	while (!eof($fastqFH)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fastqFH);
		@qualities = &$input_qual_function($qual);

		&formatSequence($fastaFH, $seq_desc, $seq);
		&formatQuality($qualFH, $seq_desc, @qualities);
	} # End of while loop

	close ($fastqFH);
	close ($fastaFH);
	close ($qualFH);
} # end of sub convertFastqToFasta

my ($fastq_file, $type, $fasta_file, $qual_file, $processors, $block);

my $result = &GetOptions("fastq|fq=s{1}" => \$fastq_file,
                         "type|t=s{1}" => \$type,
						 "fasta|f=s{1}" => \$fasta_file,
						 "qualities|q=s{1}" => \$qual_file,
						 "processors|p:i{1}" => \$processors,
						 "block|b:i{1}" => \$block);

unless ($result && defined($fastq_file) && defined($type) && defined($fasta_file) && defined($qual_file) &&
        $type =~ m/^(fastq-sanger|fastq-solexa|fastq-illumina)$/i) {
	print STDERR sprintf("\n");
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("   perl %s --fastq|-fq <fastq file> --type|-t <type> --fasta|-f <out.fas> --qualities|-q <out.qual>\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --fastq or -fq <fastq file>   : Input fastq file to be converted to standard sanger\n");
	print STDERR sprintf("                                   fasta/qualities file format\n");
	print STDERR sprintf("   --type or -t <type>           : Specify one of the fastq qualities types from below:\n");
	print STDERR sprintf("                                   'fastq-sanger'    : Sanger style fastq using PHRED and\n");
	print STDERR sprintf("                                                       ASCII offset of 33 (Range: 0 to 93)\n");
	print STDERR sprintf("                                   'fastq-solexa'    : Old solexa (prior to 1.3) quality values\n");
	print STDERR sprintf("                                                       with ASCII offset of 64 (Range: -5 to 63)\n");
	print STDERR sprintf("                                   'fastq-illumina'  : New Illumina/Solexa 1.3+ quality styles\n");
	print STDERR sprintf("                                                       with ASCII offset of 64 (Range: 0 to 63)\n");
	print STDERR sprintf("   --fasta or -f <fasta file>    : Path of output fasta file to save sequences\n");
	print STDERR sprintf("   --qualities or -q <qual file> : Path of output quality file to save quality values\n");
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

$type = lc($type);	# Convert to lower case
$processors = DEFAULT_PROCESSORS_COUNT if (!defined($processors) || $processors !~ m/^\d+$/);
$block = DEFAULT_BLOCK_SIZE if (!defined($block) || $block !~ m/^\d+$/);

my ($fastqFH, $fastaFH, $qualFH);
my $id = $$;		# Process ID
my $input_qual_function;

# Determine the fastq-variant from input
if ($type eq "fastq-sanger") {
	$input_qual_function = \&fastq_sanger_to_phred;
} # end of if statement
elsif ($type eq "fastq-solexa") {
	$input_qual_function = \&fastq_solexa_to_phred;
} # end of else if statement
elsif ($type eq "fastq-illumina") {
	$input_qual_function = \&fastq_illumina_to_phred;
} # end of else if statement
else {
	die("Invalid fastq-variant specified\n");
} # end of else statement

if ($processors == 1) {	# Uni-processor run
	&convertFastqToFasta($fastq_file, $fasta_file, $qual_file, $input_qual_function);
} # end of if statement
else {
	my %handlers;		# To keep track of different handlers
	my ($fh, $ofh);		# I/O Handllers
	my ($seq_desc, $seq, $qual_desc, $qual);

	$fh = new FileHandle();
	open ($fh, $fastq_file) or die("Cannot open input file '$fastq_file'\n");
	
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

		print $ofh &formatFastq($seq_desc, $seq, $qual_desc, $qual);
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
			&convertFastqToFasta(sprintf(".%s.%s", $p, $id), sprintf(".%s.%s.fas", $p, $id), sprintf(".%s.%s.qual", $p, $id), $input_qual_function);
			exit(0);
		} # end of else statement
	} # end of for each statement

	# Waiting for all 'forked' processes to finish
	foreach my $child (@children) {
		waitpid($child, 0);
	} # end of for each statement
	
	# Open fasta handler
	my $fastaFH = new FileHandle();
	open ($fastaFH, sprintf(">%s", $fasta_file)) or die("Cannot create final fasta file\n");
	
	# Open quality handler
	my $qualFH = new FileHandle();
	open ($qualFH, sprintf(">%s", $qual_file)) or die("Cannot create final quality file\n");

	# Merge fasta and qualities
	foreach my $p (@pids) {
		# Merge fasta
		$fh = new FileHandle();
		open ($fh, sprintf(".%s.%s.fas", $p, $id)) or die("Cannot open temporary fasta file\n");
		while (<$fh>) {
			chomp;
			print $fastaFH sprintf("%s\n", $_);
		} # end of while loop
		close ($fh);

		# Merge qualities
		$fh = new FileHandle();
		open ($fh, sprintf(".%s.%s.qual", $p, $id)) or die("Cannot open temporary quality file\n");
		while (<$fh>) {
			chomp;
			print $qualFH sprintf("%s\n", $_);
		} # end of while loop
		close ($fh);

		unlink(sprintf(".%s.%s", $p, $id));
		unlink(sprintf(".%s.%s.fas", $p, $id));
		unlink(sprintf(".%s.%s.qual", $p, $id));
	} # end of for each statement

	close ($fastaFH);
	close ($qualFH);
} # end of else statement
