#!/usr/bin/perl -w

# Author: Eddy Yeh
# 
# This script takes as command line arguments a fasta sequenes and qualities file
# and converts its corresponding values to a given fastq-variant file
#
# Change Log v0.02 - 2010.07.08 (Eddy)
#  o Changed fastq-variant quality conversion to phred scales
#
# Change Log v0.01 - 2010.06.11 (Eddy)
#  o Initial implementation of this script

use strict;
use warnings;
use FileHandle;
use POSIX qw(ceil floor);
use Getopt::Long;

use constant VERSION => "0.02 (2010.07.08)";
use constant true => 1;
use constant false => 0;

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

# Given a file handler of fasta file as argument, read and finds the next available sequence,
# returns the file handler, sequence name, and sequence nucleotides
sub nextSequence {
    my $fh = $_[0];     # file handler pointer
    my ($name, $seq);
    my $stop = false;
    my $line;

    my $pos = tell($fh);        # Obtain current file position
    while (!$stop && !eof($fh)) {
        $line = <$fh>;  # Read a line
        chomp($line);   # Remove end of line character
        $line = &trim($line);   # Remove leading and trailing white space characters

        if (length($line) != 0) {
            if ($line =~ m/^>(\S+)/) {
                if (defined($name)) {
                    seek($fh, $pos, 0); # restore position to the beginning of this line
                    $stop = true;
                } # End of if statemnet
                else {
                    $name = $1; 
                    $seq = ""; 
                } # End of else statement
            } # End of if statement
            else {
                $seq .= $line;      # Accumulate sequence nucleotides       
            } # End of else statement
        } # End of if statement

        $pos = tell($fh);   # new position before reading next line
    } # End of while loop

	if (length($seq) == 0) {	# Empty sequence
		print STDERR sprintf("\n");
		print STDERR sprintf("  ERROR: No sequence detected for '%s'\n", $name);
		print STDERR sprintf("         Program cannot proceed.\n");
		print STDERR sprintf("\n");
		exit();
	} # End of if statement

    return ($name, $seq);
} # # End of sub nextSequence

# Given a file handler of quality file as argument, read and find the next available quality,
# returns the file handler, quality name, and quality values
sub nextQuality {
    my $fh = $_[0];     # File handler pointer
    my ($name, $qual);
    my $stop = false;
    my $line;

    my $pos = tell($fh);        # Obtain current file position
    while (!$stop && !eof($fh)) {
        $line = <$fh>;  # Read a line
        chomp($line);   # Remove end of line character
        $line = &trim($line);   # Remove leading and trailing white space characters

        if (length($line) != 0) {
            if ($line =~ m/^>(\S+)/) {
                if (defined($name)) {
                    seek($fh, $pos, 0); # Restore position to the beginning of this line
                    $stop = true;
                } # End of if statement
                else {
                    $name = $1;
                    $qual = "";
                } # End of else statemenet
            } # End of if statement
            else {
                # Accumualte quality values
                if (length($qual) == 0) {
                    $qual = $line;
                } # end of if statement
                else {
                    $qual .= sprintf(" %s", $line);
                } # End of else statement
            } # End of else statement
        } # End of if statement

        $pos = tell($fh);
    } # End of while loop

    if (defined($name) && defined($qual)) { # Perform quality cleanup
        $qual = &trim($qual);       # Remove any leading and trailing white space characters
        my $index = index($qual, "  ");     # Find the first instance with 2 spaces
        while ($index >= 0) {
            $qual =~ s/  / /g;   # Replace all 2 spaces by 1 space
            $index = index($qual, "  ");
        } # End of while loop
    } # End of if statement
	
	if (length($qual) == 0) {	# Empty qualities
		print STDERR sprintf("\n");
		print STDERR sprintf("  ERROR: No quality values detected for '%s'\n", $name);
		print STDERR sprintf("         Program cannot proceed.\n");
		print STDERR sprintf("\n");
		exit();
	} # End of if statement

    return ($name, split(/\s/, $qual));
} # End of sub nextQuality

# Returns the ASCII representation of the phred quality according to its fastq-variant
sub getASCIIQuality {
	my ($type, @qual) = @_;
	
	# Convert corresponding quality to the specified format
	my $ascii_qual;
	if ($type eq "fastq-sanger") {
		$ascii_qual = &phred_to_fastq_sanger(@qual);
	} # End of if statement
	elsif ($type eq "fastq-solexa") {
		$ascii_qual = &phred_to_fastq_solexa(@qual);
	} # End of else statemnet
	else {
		$ascii_qual = &phred_to_fastq_illumina(@qual);
	} # End of else statement

	return $ascii_qual;
} # End of sub getASCIIQuality

# prints the specified fastq pieces to the file handler
sub formatFastq {
	my ($fastqFH, $seq_name, $seq, $qual_name, $qual) = @_;

	print $fastqFH sprintf("@%s\n", $seq_name);
	print $fastqFH sprintf("%s\n", $seq);
	print $fastqFH sprintf("+%s\n", $qual_name);
	print $fastqFH sprintf("%s\n", $qual);
} # End of sub formatFastq

# ---------- MAIN PROGRAM STARTS HERE ---------- #
# Command line arguments
my ($fasta_file, $quality_file, $fastq_file, $type);

my $result = &GetOptions("fasta|f=s" => \$fasta_file,
                         "qualities|q=s" => \$quality_file,
                         "fastq|fq=s" => \$fastq_file,
                         "type|t=s" => \$type);
               
unless ($result && defined($fastq_file) && defined($type) && defined($fasta_file) && defined($quality_file) &&
        $type =~ m/^(fastq-sanger|fastq-solexa|fastq-illumina)$/i) {
    print STDERR sprintf("\n");
    print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
    print STDERR sprintf("USAGE: perl %s --fasta|-f <fasta file> --qualities|-q <qual file> --type|-t <type> --fastq|-fq <fastq file>\n", $0);
	print STDERR sprintf("WHERE: --fasta or -f <fasta file>    : Input sequence file in fasta format to convert to fastq output\n");
	print STDERR sprintf("       --qualities or -q <qual file> : Input quality values of the input sequences to include in fastq output\n");
    print STDERR sprintf("       --type or -t <type>           : Specify one of the fastq qualities types from below:\n");
    print STDERR sprintf("                                       'fastq-sanger'    : Sanger style fastq using PHRED and\n");
    print STDERR sprintf("                                                           ASCII offset of 33 (Range: 0 to 93)\n");
    print STDERR sprintf("                                       'fastq-solexa'    : Old solexa (prior to 1.3) quality values\n");
    print STDERR sprintf("                                                           with ASCII offset of 64 (Range: -5 to 63)\n");
    print STDERR sprintf("                                       'fastq-illumina'  : New Illumina/Solexa 1.3+ quality styles\n");
    print STDERR sprintf("                                                           with ASCII offset of 64 (Range: 0 to 63)\n");
    print STDERR sprintf("       --fastq or -fq <fastq file>   : Path to fastq file of which output will be saved to\n");
    print STDERR sprintf("\n");
	print STDERR sprintf("VERSION: %s\n", VERSION);
	print STDERR sprintf("\n");
    exit();
} # End of unless statement

$type = lc($type);	# Conver to lower case

my ($fastaFH, $qualFH, $fastqFH);		# File handlers for input and output
my ($seq_name, $seq, $qual_name, @qual);	# sequence and quality descriptors
my (%seq_buffer, %qual_buffer);			# sequence and quality buffers in case fasta and quality file are out of order

# Open fasta handler
$fastaFH = new FileHandle();
open ($fastaFH, $fasta_file) or die("Cannot open fasta file '$fasta_file'\n");

# Open quality handler
$qualFH = new FileHandle();
open ($qualFH, $quality_file) or die("Cannot open quality file '$quality_file'\n");

# Open output handler
$fastqFH = new FileHandle();
open ($fastqFH, sprintf(">%s", $fastq_file)) or die("Cannot create fastq file '$fastq_file'\n");

# read fasta and quality handlers until either one of them is end of file
my $stop = false;
while (!$stop) {
	if (!eof($fastaFH) && !eof($qualFH)) {	# read sequence and quality and process
		($seq_name, $seq) = &nextSequence($fastaFH);
		($qual_name, @qual) = &nextQuality($qualFH);
	
		# Convert corresponding quality to the specified format
		my $ascii_qual = &getASCIIQuality($type, @qual);
		
		if ($seq_name ne $qual_name) {	# Fasta and quality file are out of order, check buffer
			if (exists $qual_buffer{$seq_name}) {	# quality for the corresponding sequence already on buffer
				&formatFastq($fastqFH, $seq_name, $seq, $seq_name, $qual_buffer{$seq_name});
				delete $qual_buffer{$seq_name};
			} # End of if statement
			else {	# Save sequence
				$seq_buffer{$seq_name} = $seq;
			} # End of else statement

			if (exists $seq_buffer{$qual_name}) {	# sequence for the corresponding quality already on buffer
				&formatFastq($fastqFH, $qual_name, $seq_buffer{$qual_name}, $qual_name, $ascii_qual);
				delete $seq_buffer{$qual_name};
			} # End of if statement
			else {	# Save quality
				$qual_buffer{$qual_name} = $ascii_qual;
			} # End of else statement
		} # End of if statement
		else {
			&formatFastq($fastqFH, $seq_name, $seq, $qual_name, $ascii_qual);
		} # End of else statement
	} # End of if statement
	elsif (!eof($fastaFH) && eof($qualFH)) {	# Only sequence present, check if its corresponding quality is in buffer
		($seq_name, $seq) = &nextSequence($fastaFH);

		if (exists $qual_buffer{$seq_name}) {	# quality for the corresponding sequence already on buffer
			&formatFastq($fastqFH, $seq_name, $seq, $seq_name, $qual_buffer{$seq_name});
			delete $qual_buffer{$seq_name};
		} # End of if statement
	} # End of else if statement
	elsif (eof($fastaFH) && !eof($qualFH)) {	# Only quality present, check if its corresponding sequence is in buffer
		($qual_name, @qual) = &nextQuality($qualFH);
		
		# Convert corresponding quality to the specified format
		my $ascii_qual = &getASCIIQuality($type, @qual);
		
		if (exists $seq_buffer{$qual_name}) {	# sequence for the corresponding quality already on buffer
			&formatFastq($fastqFH, $qual_name, $seq_buffer{$qual_name}, $qual_name, $ascii_qual);
			delete $seq_buffer{$qual_name};
		} # End of if statement
	} # End of else if statement
	else {		# fasta and quality handler reached the end
		$stop = true;
	} # End of else statement
} # End of while loop

close ($fastaFH);
close ($qualFH);
close ($fastqFH);
