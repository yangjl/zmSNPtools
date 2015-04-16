#!/usr/bin/perl -w

# Copyright (c) 2010-2011
# Cheng-Ting "Eddy" Yeh (eddyyeh@istate.edu)
# Iowa State University Research Foundation, Inc (ISURF)
# All Rights Reserved
#
# This script reads the quality values of short sequence reads from fastq variants
# (Sanger/Solexa/Illumina) and clips low quality nucleotides from the left and right
# of its sequence
#
# Change Log (v1.9) - 2011.05.24 - Eddy
#  o Small enhancement for multiprocessor trimming. --block parameter is no longer needed
#
# Change Log (v1.8) - 2011.03.04 - Eddy
#  o Added --barcode program parameter to remove 'n' number of bases at the 5' end
#    of every read thought to be the barcode sequence.
#
# Change Log (v1.7) - 2011.03.01 - Eddy
#  o Added program parameter --polyAT|-pat to enable trimming of poly A/T nucleotides
#    at the 3' end of the read. The default is to scan for >=3 consecutive A/T before
#    trimming is applied. A negative number disables this feature (i.e. do not trim
#    poly A/T tails)
#  o Replace dot (.) characters in input sequence with N's
#  o Improved/Added more options to the reported in the log file describing the quality
#    of input sequences
#  o Multiprocessor and program speed optimizations/enhancements
#
# Change Log (v1.6) - 2010.10.21 - Eddy
#  o Reimplementation of trimming reads using multi-processors. Previous implementation
#    had some problems synchronizing between processes which caused inconsitent output
#    in the final trimmed file.
#  o To overcome race conditions of trimming output, I've used temporary files to store
#    trimming results for each process which are merged at the end of the loop
#
# Change Log (v1.5.1) - 2010.08.20 - Eddy
#  o Fixed small problem in the sequence name and quality name descriptor
#
# Change Log (v1.5) - 2010.07.10 - Eddy
#  o Changed implementation to allow multiple processors for trimming if --processors|-p parameter is
#    specified. Also changed command line arguments naming for consistency purposes
#
# Change Log (v1.4) - 2010.07.08 - Eddy
#  o Changed fastq-variant quality values convertion to phred scales
#
# Change Log (v1.3) - 2010.06.24 - Eddy
#  o Changed the fastq reading implementation and avoided the use of Bioperl. Program should be much faster
#    now. This change does not affect the effects of trimming just how sequences are read from file
#
# Change Log (v1.2) - 2010.04.23 - Eddy
#  o Added trimming progress to the program
#  o Added --notrimend and --notrimmiddle options to enable/disable trimming
#    of low quality bases at the ends and middle (sliding windows) respectively.
#    Both options cannot be specified at the same time.
#
# Change Log (v1.1) - 2010.04.14 - Eddy
#  o In addition to just trim low quality bases at the ends of the sequence (front/back), I've
#    incorporated the idea of using sliding windows to detect continuous low quality
#    regions in the middle of the read, see --window option.
#    
# Change Log (v1.0) - 2010.04.12 - Eddy
#  o Initial implementation of the script
#

use strict;
use warnings;
use Getopt::Long;
use FileHandle;
use POSIX qw(ceil floor);
use Time::Local;
use Term::ANSIColor;

use constant VERSION => "1.9 (2011.05.24)";
use constant true => 1;
use constant false => 0;
use constant DEFAULT_MIN_EQ_QUALITY => 15;		# Default minimum quality threshold at the ends of the read
use constant DEFAULT_MIN_IQ_QUALITY => 15;		# Default minimum quality threshold in the middle of the read
use constant DEFAULT_MIN_LENGTH => 30;			# Default minimum length to keep after trimming
use constant DEFAULT_WINDOW_SIZE => 10;			# Default window size to use for sliding low quality regions 
use constant DEFAULT_STEP_SIZE => 1;			# Default step size per window
use constant DEFAULT_POLY_AT_SIZE => 3;			# Default poly A/T tails size before its considered
use constant DEFAULT_BARCODE_LENGTH => 0;		# Default barcode length
use constant DEFAULT_PROCESSORS_COUNT => 1;		# Default number of processors
use constant DEFAULT_QUAL_PER_LINE => 25;
use constant SLEEP_TIME => 15;					# Sleep time to wait for threads to finish execution in seconds

# Count Nucleotides
sub countNucleotides {
    my ($str, @substr) = @_; 
    my $totallength = length($str);
    
    foreach my $s (@substr) {
        $str =~ s/\Q$s//g;
    } # End of for each loop

    return ($totallength - length($str));
} # End of sub countNucleotides

# Remove comma
sub removeComma {
	my $str = $_[0];
	$str =~ s/,//g;

	return $str;
} # End of removeComma

# Formats the parametric number and adds commas every thousand, millions, etc
sub formatNumber {
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;

    return $_; 
} # End of sub formatNumber

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

# remove leading and trailing white space characters from string
sub trim {
    my $str = $_[0];
    
    if (defined($str)) {
        $str =~ s/^(\s+|\t+)//g;
        $str =~ s/(\s+|\t+)$//g;
    } # End of if statement

    return $str;
} # End of sub trim

# returns the minimum value from the parametric arguments
sub min {
	my @sorted = sort {$a <=> $b} @_;
	return shift(@sorted);
} # End of sub min

# Rounds the decimal number to integer
sub round {
    my $number = shift;
    return int($number + 0.5);
} # End of sub round

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

        if ($line =~ m/^(@\S+)/) {   # found sequence descritor
            $seq_desc = $1;
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

        if ($line =~ m/^(\+\S{0,})/) {  # found quality descriptor
            $qual_desc = $1;
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

	# Replace all dots "." w/ N
	$seq =~ tr/\./N/;

    return ($seq_desc, $seq, $qual_desc, $qual);
} # end of sub nextFastqSequence

# format the sequence so that nucleotides wrap around lines
sub formatFasta {
    my ($name, $seq) = @_; 
   	my $str;

    if ($name =~ m/^>/) {
        $str = sprintf("%s\n", $name);
    } # end of if statement
    else {
        $str = sprintf(">%s\n", $name);
    } # end of else statement

	$str .= sprintf("%s\n", $seq);

	return $str;
} # End of sub formatFasta

# format the quality values so that values wrap around lines
sub formatQuality {
    my ($name, @qual) = @_; 
	my $str;
    
    if ($name =~ m/^>/) {
        $str = sprintf("%s\n", $name);
    } # End of if statemnet
    else {
        $str = sprintf(">%s\n", $name);
    } # end of else statement

    while (scalar(@qual) > 0) {
        my @sub = splice(@qual, 0, DEFAULT_QUAL_PER_LINE);
        @sub = map { sprintf("%2s", int($_)) } @sub;
        $str .= sprintf("%s\n", join(" ", @sub));
    } # End of while loop

	return $str;
} # end of sub formatQuality

# prints the specified fastq pieces accordingly
sub formatFastq {
    my ($seq_name, $seq, $qual_name, $qual) = @_; 
	
	return sprintf("@%s\n%s\n+%s\n%s\n", $seq_name, $seq, $qual_name, $qual);
} # End of sub formatFastq

# Prints log header
sub printLogHeader {
	my ($log, $params) = @_;
	
	print $log sprintf("# %s\n", scalar(localtime(time)));
	print $log sprintf("#\n");
	print $log sprintf("# SCRIPT: %s\n", $0);
	print $log sprintf("# VERSION: %s\n", VERSION);
	print $log sprintf("#\n");
	print $log sprintf("# INPUT FASTQ: %s\n", $params->{"input"});
	print $log sprintf("# INPUT FASTQ TYPE: %s\n", $params->{"itype"});

	if ($params->{"otype"} eq "fasta") {
		print $log sprintf("# OUTPUT FASTA: %s\n", $params->{"output1"});
		print $log sprintf("# OUTPUT QUALITIES: %s\n", $params->{"output2"});
	} # End of if statemnet
	else {
		print $log sprintf("# OUTPUT FASTQ: %s\n", $params->{"output1"});
		print $log sprintf("# OUTPUT FASTQ TYPE: %s\n", $params->{"otype"});
	} # end of else statement

	if ($params->{"barcode"} == 0) {
		print $log sprintf("# BARCODE LENGTH TO REMOVE (bp): %s (no barcode removal)\n", $params->{"barcode"});
	} # end of if statement
	else {
		print $log sprintf("# BARCODE LENGTH TO REMOVE (bp): %s\n", $params->{"barcode"});
	} # end of else statement

	if ($params->{"trimend"}) {
		print $log sprintf("# MINIMUM BEGINNING/END QUALITY THRESHOLD: %s\n", $params->{"eq"});
	} # end of if statement
	else {
		print $log sprintf("# MINIMUM BEGINNING/END QUALITY THRESHOLD: DISABLED (Option: --notrimend|-note)\n");
	} # End of else statement

	if ($params->{"trimmiddle"}) {
		print $log sprintf("# WINDOW SIZE: %s\n", $params->{"window"});
		print $log sprintf("# WINDOW STEP SIZE: %s bp\n", $params->{"step"});
		print $log sprintf("# MINIMUM INTERNAL QUALITY THRESHOLD: %s\n", $params->{"iq"});
	} # end of if statemnet
	else {
		print $log sprintf("# WINDOW SIZE: DISABLED (Option: --notrimmiddle|-notm|\n");
		print $log sprintf("# WINDOW STEP SIZE: DISABLED (Option: --notrimmiddle|-notm)\n");
		print $log sprintf("# MINIMUM INTERNAL QUALITY THRESHOLD: DISABLED (Option: --notrimmiddle|-notm)\n");
	} # end of else statemnet

	if ($params->{"poly_at"} > 0) {
		print $log sprintf("# POLY A/T TRIMMING: Yes (>=%s bp)\n", $params->{"poly_at"});
	} # end of if statement
	else {
		print $log sprintf("# POLY A/T TRIMMING: No\n");
	} # end of else statemnet

	print $log sprintf("# MINIMUM READ LENGTH (POST-TRIM): %s bp\n", $params->{"minlength"});
	print $log sprintf("# PROCESSORS TO USE: %s\n", $params->{"processors"});
	print $log sprintf("#\n");
	print $log sprintf("# NOTE: Clip left and right positions are the coordinates of the read that remain\n");
	print $log sprintf("#       after trimming (good quality region). More trimming statistics at the end\n");
	print $log sprintf("#       of this file.\n");
	print $log sprintf("#\n");
	print $log sprintf("# COLUMN HEADERS: Read ID\tLength\tClip Left\tClip Right\tPost-Trim Length\n");
	print $log sprintf("\n");
} # end of sub printLogHeader

# Prints trimming statistics of the run
sub printLogStats {
	my ($log, $total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short,
        $total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT, 
		$total_bp_discarded_barcode, $minlength, $starttime) = @_;

	my $total_reads_discarded = $total_discarded_low + $total_discarded_short;
	my $total_bp_discarded = $total_bp_discarded_low + $total_bp_discarded_short + $total_bp_discarded_polyAT + $total_bp_discarded_barcode;

	print $log sprintf("\n");
	print $log sprintf("# TOTAL READS (RAW): %s\n", &formatNumber($total_reads));
	print $log sprintf("# TOTAL READS DISCARDED: %s (%2.1f%%)\n", &formatNumber($total_reads_discarded), ($total_reads_discarded / $total_reads) * 100);
	print $log sprintf("#    o READS DISCARDED DUE TO LOW QUALITY: %s (%2.1f%%)\n", &formatNumber($total_discarded_low), ($total_discarded_low / $total_reads) * 100);
    print $log sprintf("#    o READS DISCARDED DUE TO LENGTH (<%s bp): %s (%2.1f%%)\n", $minlength, 
						&formatNumber($total_discarded_short), ($total_discarded_short / $total_reads) * 100);
	print $log sprintf("# TOTAL READS REMAINING (POST-TRIM): %s (%2.1f%%)\n", &formatNumber($total_reads_post_trim), ($total_reads_post_trim / $total_reads) * 100);
	print $log sprintf("#\n");	
	print $log sprintf("# TOTAL BASE PAIRS (RAW): %s bp\n", &formatNumber($total_bp));
	print $log sprintf("# TOTAL BASE PAIRS DISCARDED: %s (%2.1f%%)\n", &formatNumber($total_bp_discarded), ($total_bp_discarded / $total_bp) * 100);
	print $log sprintf("#    o BASE PAIRS DISCARDED DUE TO BARCODE: %s (%2.1f%%)\n", &formatNumber($total_bp_discarded_barcode), ($total_bp_discarded_barcode / $total_bp) * 100);
	print $log sprintf("#    o BASE PAIRS DISCARDED DUE TO LOW QUALITY: %s (%2.1f%%)\n", &formatNumber($total_bp_discarded_low), ($total_bp_discarded_low / $total_bp) * 100);
	print $log sprintf("#    o BASE PAIRS DISCARDED DUE TO LENGTH (<%s bp): %s (%2.1f%%)\n", $minlength,
	                    &formatNumber($total_bp_discarded_short), ($total_bp_discarded_short / $total_bp) * 100);
	print $log sprintf("#    o BASE PAIRS DISCARDED DUE TO POLY A/T: %s (%2.1f%%)\n", &formatNumber($total_bp_discarded_polyAT),
	                    &formatNumber($total_bp_discarded_polyAT / $total_bp) * 100);
	print $log sprintf("# TOTAL BASE PAIRS REMAINING (POST-TRIM): %s bp (%2.2f%%)\n", &formatNumber($total_bp_post_trim), ($total_bp_post_trim / $total_bp) * 100); 
	print $log sprintf("#\n");
	print $log sprintf("# AVERAGE READ LENGTH (RAW): %d bp\n", ceil($total_bp / $total_reads));
	
	if ($total_reads_post_trim == 0) {
		print $log sprintf("# AVERAGE READ LENGTH (POST-TRIM): %d bp\n", 0);
	} # End of if statement
	else {
		print $log sprintf("# AVERAGE READ LENGTH (POST-TRIM): %d bp\n", ceil($total_bp_post_trim / $total_reads_post_trim));
	} # End of else statement

	my $endtime = timelocal(localtime(time));
	my $runtime = &formatTime($endtime - $starttime);

	print $log sprintf("#\n");
	print $log sprintf("# RUN TIME: %s\n", $runtime);
} # End of printLogStats 

# Trims input sequence and outputs in fastq format
sub trimToFastq {
	my ($params, $input, $output, $input_qual_function, $output_qual_function) = @_;
	
	my $fh = new FileHandle();
	open ($fh, $input) or die("Cannot open $input for reading\n");

	my $barcode = $params->{"barcode"};
	my $minLength = $params->{"minlength"};
	my $eq = $params->{"eq"};
	my $iq = $params->{"iq"};
	my $polyAT = $params->{"poly_at"};
	my $window = $params->{"window"};
	my $step = $params->{"step"};
	my $trimend = $params->{"trimend"};
	my $trimmiddle = $params->{"trimmiddle"};

	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s", $output)) or die("Cannot create fastq file\n");

	my ($seq_desc, $seq, $qual_desc, $qual, $clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual);

	while (!eof($fh)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fh);

		# Don't trim if sequence looks like junk
		if (&countNucleotides($seq, "N") < length($seq) / 2) {
			# Barcode checking
			$seq = substr($seq, $barcode);
			$qual = substr($qual, $barcode);

			($clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual) = 
				&trimSequence($seq, $qual, $minLength, $eq, $iq, $window, $step, $trimend, $trimmiddle, $polyAT, $input_qual_function);

			if (length($trimmed_seq) >= $minLength) {
				print $ofh &formatFastq($seq_desc, $trimmed_seq, $qual_desc, &$output_qual_function(@trimmed_qual));
			} # end of if statement
		} # end of else statemnt
	} # end of while loop
	close ($fh);
	close ($ofh);
} # End of sub trimToFastq

# Trims input sequence and outputs in fastq format with log file
sub trimToFastqLog {
	my ($params, $input, $output, $log, $input_qual_function, $output_qual_function) = @_;
	
	my $fh = new FileHandle();
	open ($fh, $input) or die("Cannot open $input for reading\n");

	my $barcode = $params->{"barcode"};
	my $minLength = $params->{"minlength"};
	my $eq = $params->{"eq"};
	my $iq = $params->{"iq"};
	my $polyAT = $params->{"poly_at"};
	my $window = $params->{"window"};
	my $step = $params->{"step"};
	my $trimend = $params->{"trimend"};
	my $trimmiddle = $params->{"trimmiddle"};

	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s", $output)) or die("Cannot create fastq file\n");

	my $lfh = new FileHandle();
	open ($lfh, sprintf(">%s", $log)) or die("Cannot create log file\n");

	&printLogHeader($lfh, $params);

	my ($seq_desc, $seq, $qual_desc, $qual, $clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual);
	my ($total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short) = (0, 0, 0, 0);
	my ($total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT) = (0, 0, 0, 0, 0);
	my ($total_bp_discarded_barcode) = 0;

	while (!eof($fh)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fh);

		# accumulate statistics
		$total_reads++;
		$total_bp += length($seq);
		
		# Don't trim if sequence looks like junk
		if (&countNucleotides($seq, "N") >= length($seq) / 2) {	# Junk
			$total_discarded_low++;
			$total_bp_discarded_low += length($seq)
		} # End of if statement
		else {
			# Barcode checking
			$seq = substr($seq, $barcode);
			$qual = substr($qual, $barcode);
			$total_bp_discarded_barcode += $barcode;
			
			($clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual) = 
				&trimSequence($seq, $qual, $minLength, $eq, $iq, $window, $step, $trimend, $trimmiddle, $polyAT, $input_qual_function);

			$total_bp_discarded_polyAT += $polya;	# Accumulate poly A/T base pairs

			if (length($trimmed_seq) >= $minLength) {
				print $ofh &formatFastq($seq_desc, $trimmed_seq, $qual_desc, &$output_qual_function(@trimmed_qual));

				# Print to log
				print $lfh sprintf("%s\t%s\t%s\t%s\t%s\n", $seq_desc, length($seq) + $barcode, $clip_left + 1 + $barcode, $clip_right + 1 + $barcode, length($trimmed_seq));

				$total_reads_post_trim++;
				$total_bp_post_trim += length($trimmed_seq);
				$total_bp_discarded_low += length($seq) - length($trimmed_seq) - $polya;
			} # end of if statement
			else {	# Discarded to to length
				$total_discarded_short++;
				$total_bp_discarded_short += length($seq) - $polya;
			} # End of else statement
		} # end of else statemnt
	} # end of while loop
	close ($fh);

	&printLogStats($lfh, $total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short,
	                $total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT, 
					$total_bp_discarded_barcode, $minLength, $params->{"start_time"});
	close ($ofh);
	close ($lfh);
} # end of sub trimToFastqLog

# Trims input sequence and outputs in fasta/quality format
sub trimToFasta {
	my ($params, $input, $output1, $output2, $input_qual_function) = @_;
	
	my $fh = new FileHandle();
	open ($fh, $input) or die("Cannot open $input for reading\n");

	my $barcode = $params->{"barcode"};
	my $minLength = $params->{"minlength"};
	my $eq = $params->{"eq"};
	my $iq = $params->{"iq"};
	my $polyAT = $params->{"poly_at"};
	my $window = $params->{"window"};
	my $step = $params->{"step"};
	my $trimend = $params->{"trimend"};
	my $trimmiddle = $params->{"trimmiddle"};

	my $ofh1 = new FileHandle();
	open ($ofh1, sprintf(">%s", $output1)) or die("Cannot create output fasta\n");

	my $ofh2 = new FileHandle();
	open ($ofh2, sprintf(">%s", $output2)) or die("Cannot create quality file\n");

	my ($seq_desc, $seq, $qual_desc, $qual, $clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual);

	while (!eof($fh)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fh);

		# Don't trim if sequence looks like junk
		if (&countNucleotides($seq, "N") < length($seq) / 2) {
			# Barcode checking
			$seq = substr($seq, $barcode);
			$qual = substr($qual, $barcode);

			($clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual) = 
				&trimSequence($seq, $qual, $minLength, $eq, $iq, $window, $step, $trimend, $trimmiddle, $polyAT, $input_qual_function);

			if (length($trimmed_seq) >= $minLength) {
				print $ofh1 &formatFasta($seq_desc, $trimmed_seq);
				print $ofh2 &formatQuality($seq_desc, @trimmed_qual);
			} # end of if statement
		} # end of else statemnt
	} # end of while loop
	close ($fh);
	close ($ofh1);
	close ($ofh2);
} # End of sub trimToFasta

# Trims input sequence and outputs in fasta/quality format with log file
sub trimToFastaLog {
	my ($params, $input, $output1, $output2, $log, $input_qual_function) = @_;
	
	my $fh = new FileHandle();
	open ($fh, $input) or die("Cannot open $input for reading\n");

	my $barcode = $params->{"barcode"};
	my $minLength = $params->{"minlength"};
	my $eq = $params->{"eq"};
	my $iq = $params->{"iq"};
	my $polyAT = $params->{"poly_at"};
	my $window = $params->{"window"};
	my $step = $params->{"step"};
	my $trimend = $params->{"trimend"};
	my $trimmiddle = $params->{"trimmiddle"};

	my $ofh1 = new FileHandle();
	open ($ofh1, sprintf(">%s", $output1)) or die("Cannot create fasta file\n");

	my $ofh2 = new FileHandle();
	open ($ofh2, sprintf(">%s", $output2)) or die("Cannot create quality file\n");

	my $lfh = new FileHandle();
	open ($lfh, sprintf(">%s", $log)) or die("Cannot create log file\n");

	&printLogHeader($lfh, $params);

	my ($seq_desc, $seq, $qual_desc, $qual, $clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual);
	my ($total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short) = (0, 0, 0, 0);
	my ($total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT) = (0, 0, 0, 0, 0);
	my ($total_bp_discarded_barcode) = 0;

	while (!eof($fh)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fh);

		# accumulate statistics
		$total_reads++;
		$total_bp += length($seq);

		# Don't trim if sequence looks like junk
		if (&countNucleotides($seq, "N") >= length($seq) / 2) {	# Junk
			$total_discarded_low++;
			$total_bp_discarded_low += length($seq)
		} # End of if statement
		else {
			# Barcode checking
			$seq = substr($seq, $barcode);
			$qual = substr($qual, $barcode);
			$total_bp_discarded_barcode += $barcode;

			($clip_left, $clip_right, $polya, $trimmed_seq, @trimmed_qual) = 
				&trimSequence($seq, $qual, $minLength, $eq, $iq, $window, $step, $trimend, $trimmiddle, $polyAT, $input_qual_function);

			$total_bp_discarded_polyAT += $polya;	# Accumulate poly A/T base pairs

			if (length($trimmed_seq) >= $minLength) {
				print $ofh1 &formatFasta($seq_desc, $trimmed_seq);
				print $ofh2 &formatQuality($seq_desc, @trimmed_qual);

				# Print to log
				print $lfh sprintf("%s\t%s\t%s\t%s\t%s\n", $seq_desc, length($seq) + $barcode, $clip_left + 1 + $barcode, $clip_right + 1 + $barcode, length($trimmed_seq));

				$total_reads_post_trim++;
				$total_bp_post_trim += length($trimmed_seq);
				$total_bp_discarded_low += length($seq) - length($trimmed_seq) - $polya;
			} # end of if statement
			else {	# Discarded to to length
				$total_discarded_short++;
				$total_bp_discarded_short += length($seq) - $polya;
			} # End of else statement
		} # end of else statemnt
	} # end of while loop
	close ($fh);

	&printLogStats($lfh, $total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short,
	                $total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT, 
					$total_bp_discarded_barcode, $minLength, $params->{"start_time"});
	close ($ofh1);
	close ($ofh2);
	close ($lfh);
} # end of sub trimToFastaLog

# Trim a sequence with the function and argument passed
sub trimSequence {
	my ($sequence, $quality, $minLength, $eq, $iq, $window, $step, $trimend, $trimmiddle, $polyAT, $input_qual_function) = @_;

	my @seq = split(//, $sequence);
	my @qual = &$input_qual_function($quality);
			
	my $clip_left = 0;
	my $clip_right = scalar(@qual) - 1;

	if ($trimend) {	# Trim ends if enabled
		# Trim low quality nucleotides < End Quality Threshold  at the beginning of the sequence
		# Stop at the first nucleotide that is encontered with quality >= End Qualilty Threshold
		while ($qual[$clip_left] < $eq && $clip_left < $clip_right) {
			$clip_left++;
		} # End of while loop

		# Trim low quality nucleotides < End Quality Threshold at the end of the sequence
		# Stop at the first nucleotide that is encontered with quality >= End Quality Threshold
		while ($qual[$clip_right] < $eq && $clip_right > $clip_left) {
			$clip_right--;
		} # End of while loop
	} # End of if statement

	if ($trimmiddle) {	# Trim middle by sliding windows if enabled
		# For the remaining middle part of the sequence, scan every window stepping each 'step'
		# size to determine the average quality value. If average quality value is < Internal Quality Threshold,
		# trim the remaining of the sequence
		my $continue = true;
		for (my $i=$clip_left + 1; $i < $clip_right && $continue; $i = $i + $step) {
			my $sum = 0;
			my $times = 0;
			for (my $j = 0; $j < $window && $j + $i < $clip_right; $j++) {
				$sum += $qual[$j + $i];
				$times++;
			} # End of for loop

			if ($times == 0) {	# Nothing was trimmed
				$continue = false;
			} # End of if statement
			elsif (&round($sum / $times) < $iq) { 	# first low quality window
				$continue = false;
				# update clip right coordinate
				$clip_right = $i - 1;
			} # End of if statement
		} # End of for loop
	} # End of if statement

	# Extract clipped sequences and qualities
	@seq = splice(@seq, $clip_left, $clip_right - $clip_left + 1);
	@qual = splice(@qual, $clip_left, $clip_right - $clip_left + 1);

	$sequence = join("", @seq);

	# PolyA/T Trimming
	my $polya = 0;		# Holds the total number of nucleotides trimmed due to poly A/T
	if ($polyAT >= 0) {
		if ($sequence =~ m/([A|T]{\Q$polyAT\E,})$/i) {
			$polya = length($1);
			
			# Update sequence
			$sequence = substr($sequence, 0, length($sequence) - $polya);
			@qual = splice(@qual, 0, length($sequence));
			$clip_right -= $polya;
		} # End of if statement
	} # end of if statement

	return ($clip_left, $clip_right, $polya, $sequence, @qual);
} # end of trimSequence

# I/O parameters
my ($input, $itype, @output, $otype, $logfile, $help, $processors);
my %params;
my $total_reads = 0;
my $total_bp = 0;
my $total_reads_post_trim = 0;
my $total_bp_post_trim = 0;

# ------------------ MAIN PROGRAM STARTS HERE -------------------- #
my $result = &GetOptions("input|i=s" => \$input,
                         "inputtype|itype|it=s" => \$itype,
						 "output|o=s{1,2}" => \@output,
						 "outputtype|otype|ot=s" => \$otype,
						 "barcode:i{1}" => \$params{"barcode"},
						 "endquality|endqual|eq:i" => \$params{"eq"},
						 "minimum|m:i" => \$params{"minlength"},
						 "window|w:i{1}" => \$params{"window"},
						 "internalquality|intqual|iq:i" => \$params{"iq"},
						 "step|s:i" => \$params{"step"},
						 "trimend|te!" => \$params{"trimend"},
						 "trimmiddle|tm!" => \$params{"trimmiddle"},
						 "polyAT|pat:i{1}" => \$params{"poly_at"},
						 "processors|proc|p:i{1}" => \$processors,
						 "log|l:s" => \$logfile,
						 "help|h" => \$help);

unless ($result && defined($input) && defined($itype) && defined($otype) && (scalar(@output) >= 1 && scalar(@output) <= 2) && 
        $itype =~ m/^(fastq-sanger|fastq-solexa|fastq-illumina)$/i && $otype =~ m/^(fasta|fastq-sanger|fastq-solexa|fastq-illumina)$/i) {
	print STDERR sprintf("\n");
	print STDERR colored("*** INCORRECT NUMBER OF ARGUMENTS ***", "bold red ON_BLACK") . "\n" if (!defined($help));
	print STDERR colored("USAGE:", "bold yellow ON_BLACK") . "\n";
	print STDERR colored(sprintf("   perl %s --input|-i <fastq file> --inputtype|--itype|-it <type> [OPTIONS]", $0), "bold green ON_BLACK") . "\n";
	print STDERR sprintf("\n");
	print STDERR colored("WHERE:", "bold yellow ON_BLACK") . "\n";
	print STDERR sprintf("   --input|-i <fastq file>                  : Input fastq file to be converted to standard sanger\n");
	print STDERR sprintf("                                              fasta/qualities file format\n");
	print STDERR sprintf("   --inputtype|--itype|-it <type>           : Specify one of the fastq qualities types from below:\n");
	print STDERR sprintf("                                              'fastq-sanger'    : Sanger style fastq with ASCII encoding from 33 to 105\n");
	print STDERR sprintf("                                                                  (ASCII-33) ; PHRED quality score from 0 to 93\n");
	print STDERR sprintf("                                              'fastq-solexa'    : Solexa quality w/ Cassava pipeline prior to 1.3\n");
	print STDERR sprintf("                                                                  with ASCII encoding from 59 to 104 (ASCII-64); PHRED\n");
	print STDERR sprintf("                                                                  quality score from -5 to 40\n");
	print STDERR sprintf("                                              'fastq-illumina'  : Illumina/Solexa Cassava pipeline 1.3+ with ASCII encoding\n");
	print STDERR sprintf("                                                                  from 64 to 104 (ASCII-64) ; PHRED quality score from 0 to 40\n");
	print STDERR sprintf("\n");
	print STDERR colored("TRIMMING OPTIONS:", "bold yellow ON_BLACK") . "\n";
	print STDERR sprintf("   --barcode <barcode length>               : Specify the number of bases to remove at the 5' end thought to be\n");
	print STDERR sprintf("                                              sequence barcodes [DEFAULT: %s]\n", DEFAULT_BARCODE_LENGTH);
	print STDERR sprintf("   --endquality|--endqual|-eq <min. qual>   : Minimum quality threshold to consider for trimming in\n");
	print STDERR sprintf("                                              decimal format (Range: 0 to 40) at the begining and ends\n");
	print STDERR sprintf("                                              of the read [DEFAULT: %s]\n", DEFAULT_MIN_EQ_QUALITY);
	print STDERR sprintf("   --minimum|-m <min. length>               : Minimum length of sequence after trimming. [DEFAULT: %s]\n", DEFAULT_MIN_LENGTH);
	print STDERR sprintf("   --internalquality|intqual|iq <min. qual> : Minimum average quality threshold to consider for trimming in\n");
	print STDERR sprintf("                                              decimal format (Range: 0 to 40) per window [DEFAULT: %s]\n", DEFAULT_MIN_IQ_QUALITY);
	print STDERR sprintf("   --window|-w <size>                       : Specify the number of bases to use as the sliding window size\n");
	print STDERR sprintf("                                              in the middle of the read [DEFAULT: %s]\n", DEFAULT_WINDOW_SIZE);
	print STDERR sprintf("   --step|-s <step size>                    : Specify the number of bases to the step per window [DEFAULT: %s]\n", DEFAULT_STEP_SIZE);
	print STDERR sprintf("   --trimend|-te|--notrimend|-note          : Enables (--trimend or -te) or disables (--notrimend or -note) the\n");
	print STDERR sprintf("                                              trimming of low quality nucleotides at the beginning and end of\n");
	print STDERR sprintf("                                              of the read if quality value is below --endquality [DEFAULT: --trimend]\n");
	print STDERR sprintf("   --trimmiddle|-tm|--notrimmiddle|-notm    : Enables (--trimmiddle or -tm) or disables (--notrimmiddle or -notm) the\n");
	print STDERR sprintf("                                              trimming of low quality nucleotides detected by sliding windows (--window\n");
	print STDERR sprintf("                                              and --step) across the entire read [DEFAULT: --trimmiddle]\n");
	print STDERR sprintf("   --polyAT|-pat <length>                   : Specify the number of consecutive poly A/T nucleotides allowed before is\n");
	print STDERR sprintf("                                              considered as polyadenylation tail. A negative number disables trimming of\n");
	print STDERR sprintf("                                              poly A/T tails [DEFAULT: %s]\n", DEFAULT_POLY_AT_SIZE);
	print STDERR sprintf("   --processors|--proc|-p <num processors>  : If multiple processors is available, specify the number of \"virtual\"\n");
	print STDERR sprintf("                                              processes to use for trimming. Essentially, this option divides the input\n");
	print STDERR sprintf("                                              sequences into batches which are processed independently and merging them\n");
	print STDERR sprintf("                                              to form the final output [DEFAULT: %s]\n", DEFAULT_PROCESSORS_COUNT);
	print STDERR sprintf("\n");
	print STDERR colored("OUTPUT OPTIONS:", "bold yellow ON_BLACK") . "\n";
	print STDERR sprintf("   --outputtype|-otype|-ot <type>           : Specify one of the following output formats from below:\n");
	print STDERR sprintf("                                              'fasta'           : Outputs the trimmed sequences in fasta/quality\n");
	print STDERR sprintf("                                                                  formats. If this option is specified, --output|-o\n");
	print STDERR sprintf("                                                                  must include TWO (2) file names to save trimmed\n");
	print STDERR sprintf("                                                                  sequences and qualities respectively\n");
	print STDERR sprintf("                                              'fastq-sanger'    : Outputs the trimmed sequences and quality values\n");
	print STDERR sprintf("                                                                  in Sanger style fastq format with ASCII offset of\n");
	print STDERR sprintf("                                                                  33. If this option is specified, --output|-o must\n");
	print STDERR sprintf("                                                                  include ONE (1) file name to save trimmed sequences\n");
	print STDERR sprintf("                                                                  and qualities\n");
	print STDERR sprintf("                                              'fastq-solexa'    : Outputs the trimmed sequences and quality values\n");
	print STDERR sprintf("                                                                  in Solexa Cassava pipeline prior to 1.3 fastq format with\n");
	print STDERR sprintf("                                                                  ASCII offset of 64. If this option is specified, --output|-o\n");
	print STDERR sprintf("                                                                  must include ONE (1) file name to save trimmed sequences\n");
	print STDERR sprintf("                                                                   and qualities\n");
	print STDERR sprintf("                                               'fastq-illumina' : Outputs the trimmed sequences and quality values\n");
	print STDERR sprintf("                                                                  in fastq Illumina format (1.3+) with ASCII offset of 64.\n");
	print STDERR sprintf("                                                                  If this option is specified, --output|-o must include\n");
	print STDERR sprintf("                                                                  ONE (1) file name to save trimmed sequences and qualities\n");
	print STDERR sprintf("   --output|-o <file 1> [<file 2>]          : File names to save trimmed sequences. If --otype is 'fasta', specify TWO (2)\n");
	print STDERR sprintf("                                              files, one for trimmed sequences and the second one for trimmed qualities.\n");
	print STDERR sprintf("                                              Specify only ONE (1) file for all other output formats\n");
	print STDERR sprintf("   --log|-l <log file>                      : If this option is specified, output trim information to the file specified\n");
	print STDERR sprintf("\n");
	print STDERR colored("VERSION: ", "bold yellow ON_BLACK") . sprintf("%s\n", VERSION);
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

# Save input type format
$itype = lc($itype);	# Convert to lower case
$otype = lc($otype);	# Convert to lower case

# Assigning default values
$params{"barcode"} = DEFAULT_BARCODE_LENGTH if (!defined($params{"barcode"}) || !exists $params{"barcode"} ||  $params{"barcode"} < 0);
$params{"eq"} = DEFAULT_MIN_EQ_QUALITY if (!defined($params{"eq"}) || !exists $params{"eq"} ||  $params{"eq"} < 0);
$params{"minlength"} = DEFAULT_MIN_LENGTH if (!defined($params{"minlength"}) || !exists $params{"minlength"} || $params{"minlength"} < 0);
$params{"iq"} = DEFAULT_MIN_IQ_QUALITY if (!defined($params{"iq"}) || !exists $params{"iq"} || $params{"iq"} < 0);
$params{"window"} = DEFAULT_WINDOW_SIZE if (!defined($params{"window"}) || !exists $params{"window"} || $params{"window"} < 0);
$params{"step"} = DEFAULT_STEP_SIZE if (!defined($params{"step"}) || !exists $params{"step"} || $params{"step"} < 0);
$params{"poly_at"} = DEFAULT_POLY_AT_SIZE if (!defined($params{"poly_at"}) || !exists $params{"poly_at"});
$params{"trimend"} = true if (!defined($params{"trimend"}) || !exists $params{"trimend"});
$params{"trimmiddle"} = true if (!defined($params{"trimmiddle"}) || !exists $params{"trimmiddle"});
$processors = DEFAULT_PROCESSORS_COUNT if (!defined($processors) || $processors !~ m/^(\d+)$/);

# Saving more values
$params{"input"} = $input;
$params{"itype"} = $itype;
if ($otype eq "fasta") {
	$params{"output1"} = $output[0];
	$params{"output2"} = $output[1];
} # end of if statement
else {
	$params{"output1"} = $output[0];
} # end of else statement
$params{"otype"} = $otype;
$params{"processors"} = $processors;

# Start time of program
$params{"start_time"} = timelocal(localtime(time));

if ($otype eq "fasta" && scalar(@output) != 2) {
	print STDERR sprintf("\n");
	print STDERR sprintf("You must specify TWO (2) file names in --output|-o option to save trimmed sequences and qualities\n");
	print STDERR sprintf("respectively if --otype|-ot is 'fasta'\n");
	print STDERR sprintf("\n");
	exit();
} # End of if statement
elsif ($otype ne "fasta" && scalar(@output) != 1) {
	print STDERR sprintf("\n");
	print STDERR sprintf("You must specify ONE (1) file name in --output|-o option to save trimmed sequences/qualities in fastq\n");
	print STDERR sprintf("format if --otype|-ot is either 'fastq-sanger', 'fastq-solexa', or 'fastq-illumina'\n");
	print STDERR sprintf("\n");
	exit();
} # End of else statement
elsif ($params{"trimend"} == false && $params{"trimmiddle"} == false) {
	print STDERR sprintf("\n");
	print STDERR sprintf("You must specify either to trim ends, trim middle, or both. --notrimend and --notrimmiddle cannot be\n");
	print STDERR sprintf("both specified since it defeats the purpose of trimming.\n");
	print STDERR sprintf("\n");
	exit();
} # End of else if statement

if (defined($logfile) && length($logfile) == 0) {
	print STDERR sprintf("\n");
	print STDERR sprintf("You must specify a valid file path for --log|-l option\n");
	print STDERR sprintf("\n");
	exit();
} # end of if statement

# -------------------- MAIN PROGRAM STARTS HERE -------------------- #
my $fastqFH;								# handler to read fastq sequences
my $ofh;									# output handler
my $id = $$;								# job id
my ($seq_desc, $seq, $qual_desc, $qual);	# I/O variables

# Determine input and output quality functions
my $input_qual_function;
if ($itype =~ m/^fastq-sanger$/i) {
	$input_qual_function = \&fastq_sanger_to_phred;
} # end of if statement
elsif ($itype =~ m/^fastq-solexa$/i) {
	$input_qual_function = \&fastq_solexa_to_phred;
} # end of else statement
else {
	$input_qual_function = \&fastq_illumina_to_phred;
} # end of else statement

my $output_qual_function;
if ($otype =~ m/^fastq-/i) {	# Fastq output
	if ($otype =~ m/^fastq-sanger$/i) {
		$output_qual_function = \&phred_to_fastq_sanger;
	} # end of if statement
	elsif ($otype =~ m/^fastq-solexa$/i) {
		$output_qual_function = \&phred_to_fastq_solexa;
	} # end of elsif statement
	elsif ($otype =~ m/^fastq-illumina$/i) {
		$output_qual_function = \&phred_to_fastq_illumina;
	} # end of else statement
	else {
		die("This shouldn't happen, invalid otype argument\n");
	} # End of else statement
} # end of if statement

if ($processors == 1) {	# Using single processor
	if ($otype =~ m/^fastq-/i) {	# Fastq output
		if (!defined($logfile)) {
			&trimToFastq(\%params, $input, $output[0], $input_qual_function, $output_qual_function);
		} # end of if statement
		else {
			&trimToFastqLog(\%params, $input, $output[0], $logfile, $input_qual_function, $output_qual_function);
		} # end of else statement
	} # end of if statement
	else {	# Fasta/Qual output
		if (!defined($logfile)) {
			&trimToFasta(\%params, $input, $output[0], $output[1], $input_qual_function);
		} # end of if statement
		else {
			&trimToFastaLog(\%params, $input, $output[0], $output[1], $logfile, $input_qual_function);
		} # end of else statemetn
	} # end of else statmenet
} # end of if statement
else {	# Multi-processor triming, split input file and virtually parallelize the job
	my ($acc_reads, $acc_bp) = (0, 0);	# Accumulators
	$fastqFH = new FileHandle();
	open ($fastqFH, $input) or die("cannot input file $input\n");

	my %handlers;	# To keep track of different handlers for each splitted file
	my $total_sequences = `LANG=C wc -l $input`;

	if ($total_sequences =~ m/^(\d+)\s+\S+/) {
		$total_sequences = $1 / 4;
	} # end of if statement
	else {
		die("Cannot determine the total number of sequences in file $input\n");
	} # end of else statement

	my $perFile = ceil($total_sequences / $processors);

	my $pnum = 0;
	my $count = 0;
	while (!eof($fastqFH)) {
		($seq_desc, $seq, $qual_desc, $qual) = &nextFastqSequence($fastqFH);

		if ($count == 0 || $count % $perFile == 0) {
			close ($ofh) if (defined($ofh));
			$ofh = new FileHandle();
			open ($ofh, sprintf(">.%s.%s", $pnum, $id)) or die("cannot create temporary file\n");
			$handlers{$pnum} = $ofh;
			$pnum++;
		} # end of if statemnet

		print $ofh &formatFastq($seq_desc, $seq, $qual_desc, $qual);
		$count++;

		# Accumulate
		$acc_reads++;
		$acc_bp += length($seq);
	} # end of while loop
	close ($fastqFH);
	close ($ofh) if (defined($ofh));
	
	my @pids = sort {$a <=> $b} keys %handlers;

	# Close all open file handlers
	foreach my $p (@pids) {
		close($handlers{$p});
	} # end of for each statement

	# Initializing multiprocessor run
	my @children;
	foreach my $p (@pids) {
		my $pid = fork();		# Spawn a new process

		if ($pid) {	# Parent process
			push(@children, $pid);
		} # end of if statement
		else {	# Child process performs trimming
			if ($otype =~ m/^fastq-/i) {
				if (!defined($logfile)) {
					&trimToFastq(\%params, sprintf(".%s.%s", $p, $id), sprintf(".%s.%s.seq", $p, $id), $input_qual_function, $output_qual_function);
				} # End of if statement
				else {
					&trimToFastqLog(\%params, sprintf(".%s.%s", $p, $id), sprintf(".%s.%s.seq", $p, $id), sprintf(".%s.%s.log", $p, $id), $input_qual_function, $output_qual_function);
				} # end of else statemnet
			} # End of if statement
			else {
				if (!defined($logfile)) {
					&trimToFasta(\%params, sprintf(".%s.%s", $p, $id), sprintf(".%s.%s.seq", $p, $id), sprintf(".%s.%s.qual", $p, $id), $input_qual_function);
				} # End of if statement
				else {
					&trimToFastaLog(\%params, sprintf(".%s.%s", $p, $id), sprintf(".%s.%s.seq", $p, $id), sprintf(".%s.%s.qual", $p, $id), sprintf(".%s.%s.log", $p, $id), $input_qual_function);
				} # end of else statement
			} # End of else statement
		
			exit(0);	# Exit this child process
		} # end of else statement
	} # end of for each statement

	# Waiting for all 'forked' processes to finish
	foreach my $child (@children) {
		waitpid($child, 0);
	} # End of for each statement

	# Merge output from all processes
	my $ofh1 = new FileHandle();
	open ($ofh1, sprintf(">%s", $output[0])) or die("Cannot create final output file for sequences\n");

	my $ofh2;
	if ($otype eq "fasta") {
		$ofh2 = new FileHandle();
		open ($ofh2, sprintf(">%s", $output[1])) or die("Cannot create final output file for qualities\n");
	} # end of if statement

	my $ofh3;
	my ($total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short) = (0, 0, 0, 0, 0);
	my ($total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT) = (0, 0, 0, 0, 0);
	my ($total_bp_discarded_barcode) = 0;
	my ($acc_reads_post_trim, $acc_bp_post_trim) = (0, 0);
	if (defined($logfile)) {
		$ofh3 = new FileHandle();
		open ($ofh3, sprintf(">%s", $logfile)) or die("Cannot create final log file\n");

		&printLogHeader($ofh3, \%params);
	} # end of if statement

	# Merging sequence file
	foreach my $p (@pids) {
		$fastqFH = new FileHandle();
		open ($fastqFH, sprintf(".%s.%s.seq", $p, $id)) or die("Cannot open temporary file\n");
		while (<$fastqFH>) {
			chomp;
			print $ofh1 sprintf("%s\n", $_);
		} # End of while loop
		close ($fastqFH);
	
		if ($otype eq "fasta") {	# Need to merge qualities
			$fastqFH = new FileHandle();
			open ($fastqFH, sprintf(".%s.%s.qual", $p, $id)) or die("Cannot open temporary file\n");
			while (<$fastqFH>) {
				chomp;
				print $ofh2 sprintf("%s\n", $_);
			} # End of while loop
			close ($fastqFH);
		} # End of if statement

		if (defined($logfile)) {	# Need to merge log files
			$fastqFH = new FileHandle();
			open ($fastqFH, sprintf(".%s.%s.log", $p, $id)) or die("Cannot open temporary file\n");
			while (<$fastqFH>) {
				chomp;
				if (length($_) != 0) {
					if ($_ =~ m/^# TOTAL READS \(RAW\): (\S+)/) {
						$total_reads += &removeComma($1);	
					} # end of if statement
					elsif ($_ =~ m/^# TOTAL READS REMAINING \(POST-TRIM\): (\S+)/) {
						$total_reads_post_trim += &removeComma($1);
					} # end of else if statement
					elsif ($_ =~ m/^#    o READS DISCARDED DUE TO LOW QUALITY: (\S+)/) {
						$total_discarded_low += &removeComma($1);
					} # end of elsif statement
					elsif ($_ =~ m/^#    o READS DISCARDED DUE TO LENGTH \(<\Q$params{"minlength"}\E bp\): (\S+)/) {
						$total_discarded_short += &removeComma($1);
					} # End of elsif statement
					elsif ($_ =~ m/^# TOTAL READS REMAINING \(POST-TRIM\): (\S+)/) {
						$total_reads_post_trim += &removeComma($1);
					} # End of else if statement
					elsif ($_ =~ m/^# TOTAL BASE PAIRS \(RAW\): (\S+)/) {
						$total_bp += &removeComma($1);
					} # end of else if statement
					elsif ($_ =~ m/^# TOTAL BASE PAIRS REMAINING \(POST-TRIM\): (\S+)/) {
						$total_bp_post_trim += &removeComma($1);
					} # End of if statement
					elsif ($_ =~ m/^#    o BASE PAIRS DISCARDED DUE TO BARCODE: (\S+)/) {
						$total_bp_discarded_barcode += &removeComma($1);
					} # end of else if statement
					elsif ($_ =~ m/^#    o BASE PAIRS DISCARDED DUE TO LOW QUALITY: (\S+)/) {
						$total_bp_discarded_low += &removeComma($1);
					} # end of else if statement
					elsif ($_ =~ m/^#    o BASE PAIRS DISCARDED DUE TO LENGTH \(<\Q$params{"minlength"}\E bp\): (\S+)/) {
						$total_bp_discarded_short += &removeComma($1);
					} # end of else if statement
					elsif ($_ =~ m/^#    o BASE PAIRS DISCARDED DUE TO POLY A\/T: (\S+)/) {
						$total_bp_discarded_polyAT += &removeComma($1);
					} # end of else if statement
					elsif ($_ !~ m/^#/) {
						my @fields = split(/\t/, $_);
						$acc_reads_post_trim++;
						$acc_bp_post_trim += $fields[$#fields];
						print $ofh3 sprintf("%s\n", $_);
					} # end of if statement
				} # End of if statement
			} # End of while loop
		} # end of if statement

		unlink(sprintf(".%s.%s", $p, $id));
		unlink(sprintf(".%s.%s.seq", $p, $id));
		unlink(sprintf(".%s.%s.qual", $p, $id)) if ($otype eq "fasta");
		unlink(sprintf(".%s.%s.log", $p, $id)) if (defined($logfile));
	} # end of for each statement

	close ($ofh1);
	close ($ofh2) if (defined($ofh2));

	if (defined($ofh3)) {
		# Display warning if numbers don't match
		if ($acc_reads != $total_reads || $acc_bp != $total_bp || $acc_reads_post_trim != $total_reads_post_trim || $acc_bp_post_trim != $total_bp_post_trim) {
			print STDERR sprintf("WARNING\n");
			print STDERR sprintf("   o Some processes in this run did not complete successfully. Some sequences might be missing from the\n");
			print STDERR sprintf("     final trimming output\n");
		} # End of if statement

		&printLogStats($ofh3, $total_reads, $total_reads_post_trim, $total_discarded_low, $total_discarded_short,
        			   $total_bp, $total_bp_post_trim, $total_bp_discarded_low, $total_bp_discarded_short, $total_bp_discarded_polyAT, 
					   $total_bp_discarded_barcode, $params{"minlength"}, $params{"start_time"});
		close ($ofh3);
	} # end of if statement
} # end of else statement
# -------------------- END OF MAIN PROGRAM ------------------ #
