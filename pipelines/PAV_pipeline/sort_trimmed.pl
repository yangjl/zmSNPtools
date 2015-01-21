#!/usr/bin/perl -w

# Author: Eddy Yeh
# Schnable Laboratory, Iowa State University

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil floor);
use FileHandle;

use constant VERSION => "0.02 (2013.03.08)";
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

	my @tmp = split(/\s/, $seq_desc);
	$seq_desc = shift(@tmp);

	if (length($qual_desc) > 1) {
		@tmp = split(/\s/, $qual_desc);
		$qual_desc = shift(@tmp);
	} # end of if statement

    return ($seq_desc, $seq, $qual_desc, $qual);
} # end of sub nextFastqSequence

sub printFastqSequence {
	my ($outFH, $seqdesc, $seq, $qualdesc, $qual, $mate) = @_;

	if ($seqdesc =~ m/^@/) {
		if ($seqdesc !~ m/\/(1|2)$/) {
			print $outFH sprintf("%s/%s\n", $seqdesc, $mate);
		} # end of if statement
		else {
			print $outFH sprintf("%s\n", $seqdesc);
		} # End of else statement
	} # End of if statemnet
	else {
		if ($seqdesc !~ m/\/(1|2)$/) {
			print $outFH sprintf("@%s/%s\n", $seqdesc, $mate);
		} # end of if statement
		else {
			print $outFH sprintf("@%s\n", $seqdesc);
		} # End of else statement
	} # End of else statemnt
	
	print $outFH sprintf("%s\n", $seq);

	if ($qualdesc =~ m/^\+/) {
		if (length($qualdesc) > 1 && $qualdesc !~ m/\/(1|2)$/) {
			print $outFH sprintf("%s/%s\n", $qualdesc, $mate);
		} # end of if statement
		else {
			print $outFH sprintf("%s\n", $qualdesc);
		} # end of else statmenet
	} # End of if statemnet
	else {
		if (length($qualdesc) > 1 && $qualdesc !~ m/\/(1|2)$/) {
			print $outFH sprintf("+%s/%s\n", $qualdesc, $mate);
		} # End of if statemnet
		else {
			print $outFH sprintf("+%s\n", $qualdesc);
		} # end of else statement
	} # End of else statement

	print $outFH sprintf("%s\n", $qual);
} # end of printFastqSequence

my ($file1, $file2, $prefix, $singletons);

my $result = &GetOptions("file1|f1=s{1}" => \$file1,
                         "file2|f2=s{1}" => \$file2,
						 "prefix|p=s{1}" => \$prefix,
						 "singletons|s!" => \$singletons);

unless ($result && defined($file1) && defined($file2) && defined($prefix)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("   perl %s --file1|-f1 <file1.fq> --file2|-f2 <file2.fq> --prefix <out> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --file1|-f1 <file1.fq>              : Path to read #1 file in fastq format\n");
	print STDERR sprintf("   --file2|-f2 <file2.fq>              : Path to read #2 file in fastq format\n");
	print STDERR sprintf("   --prefix <out>                      : Prefix string for output files\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("   --singletons|-s|--nosingletons|-nos : Separate singleton reads as independent files\n");
	print STDERR sprintf("                                         Default behaviour is to separate singleton reads (--singletons).\n");
	print STDERR sprintf("                                         If either one of the reads from #1 or #2 is not present\n");
	print STDERR sprintf("                                         an empty sequence will be placed as their corresponding read\n");
	print STDERR sprintf("                                         otherwise if --singletons is specified, corresponding\n");
	print STDERR sprintf("                                         <prefix>-singletons-1.fq and <prefix>-singletons-2.fq files will\n");
	print STDERR sprintf("                                         be generated\n");
	print STDERR sprintf("\n");
	exit();
} # end of else statement

# Assign default singletons behavior
$singletons = true if (!defined($singletons));

my (%buffer1, %buffer2);
my ($fh1, $fh2);
my ($out1, $out2);

# Open file handlers for reading
open ($fh1, $file1) or die("Cannot open file 1 for reading\n");
open ($fh2, $file2) or die("Cannot open file 2 for reading\n");

# Open file handlers for writing
open ($out1, sprintf(">%s-paired-1.fq", $prefix)) or die("Cannot create output file 1\n");
open ($out2, sprintf(">%s-paired-2.fq", $prefix)) or die("Cannot create outptu file 2\n");

my ($seqdesc1, $seq1, $qualdesc1, $qual1);
my ($seqdesc2, $seq2, $qualdesc2, $qual2);
my ($name1, $name2);

my $continue = true;
while ($continue) {
	if (!eof($fh1)) {
		($seqdesc1, $seq1, $qualdesc1, $qual1) = &nextFastqSequence($fh1);
		
		if ($seqdesc1 =~ m/\/(1|2)$/) {
			$name1 = $`;
		} # End of if statemnet
		else {
			$name1 = $seqdesc1;
		} # End of else statement
		
		if (length($seq1) == 0) {
			undef($name1);
		} # End of if statement
	} # end of if statemnet
	else {
		undef($name1);
	} # end of if statement

	if (!eof($fh2)) {
		($seqdesc2, $seq2, $qualdesc2, $qual2) = &nextFastqSequence($fh2);
	
		if ($seqdesc2 =~ m/\/(1|2)$/) {
			$name2 = $`;
		} # End of if statemnet
		else {
			$name2 = $seqdesc2;
		} # End of else statemnet

		if (length($seq2) == 0) {
			undef($name2);
		} # End of if statement
	} # End of if statemnet
	else {
		undef($name2);
	} # End of else statement

	# Checking
	if (defined($name1) && defined($name2) && $name1 eq $name2) {	# Found the corresponding mate
		&printFastqSequence($out1, $seqdesc1, $seq1, $qualdesc1, $qual1, 1);
		&printFastqSequence($out2, $seqdesc2, $seq2, $qualdesc2, $qual2, 2);
	} # end of if statemnet
	else {
		if (defined($name1)) {
			if (exists $buffer2{$name1}) {	# mate of file handler 1 was already read before
				&printFastqSequence($out1, $seqdesc1, $seq1, $qualdesc1, $qual1, 1);
				&printFastqSequence($out2, $buffer2{$name1}->{"seqdesc"}, $buffer2{$name1}->{"seq"},
				                           $buffer2{$name1}->{"qualdesc"}, $buffer2{$name1}->{"qual"}, 2);
				delete $buffer2{$name1};	# no longer needed
			} # End of if statement
			else {	# save this entry on buffer 1
				$buffer1{$name1}->{"seqdesc"} = $seqdesc1;
				$buffer1{$name1}->{"seq"} = $seq1;
				$buffer1{$name1}->{"qualdesc"} = $qualdesc1;
				$buffer1{$name1}->{"qual"} = $qual1;
			} # End of else statemnet
		} # End of if statemnet

		if (defined($name2)) {
			if (exists $buffer1{$name2}) {	# mate of file handler 2 was already read before
				&printFastqSequence($out1, $buffer1{$name2}->{"seqdesc"}, $buffer1{$name2}->{"seq"},
				                           $buffer1{$name2}->{"qualdesc"}, $buffer1{$name2}->{"qual"}, 1);
				&printFastqSequence($out2, $seqdesc2, $seq2, $qualdesc2, $qual2, 2);
				delete $buffer1{$name2};	# no longer needed
			} # End of if statement
			else {	# save this entry on buffer 2
				$buffer2{$name2}->{"seqdesc"} = $seqdesc2;
				$buffer2{$name2}->{"seq"} = $seq2;
				$buffer2{$name2}->{"qualdesc"} = $qualdesc2;
				$buffer2{$name2}->{"qual"} = $qual2;
			} # End of else statement
		} # End of if statement
	} # End of else statement

	$continue = false if (eof($fh1) && eof($fh2));
} # End of while loop
close ($fh1);
close ($fh2);

if ($singletons) {
	close ($out1);
	close ($out2);

	# print singletons
	my @keys = keys %buffer1;
	if (scalar(@keys) > 0) {
		open ($out1, sprintf(">%s-singletons-1.fq", $prefix)) or die("Cannot create singletons file 1\n");
		foreach my $k (@keys) {
			&printFastqSequence($out1, $buffer1{$k}->{"seqdesc"}, $buffer1{$k}->{"seq"},
									   $buffer1{$k}->{"qualdesc"}, $buffer1{$k}->{"qual"}, 1);
		} # End of for each statement
		close ($out1);
	} # End of if statemnet

	@keys = keys %buffer2;
	if (scalar(@keys) > 0) {
		open ($out2, sprintf(">%s-singletons-2.fq", $prefix)) or die("Cannot create singletons file 2\n");
		foreach my $k (@keys) {
			&printFastqSequence($out2, $buffer2{$k}->{"seqdesc"}, $buffer2{$k}->{"seq"},
									   $buffer2{$k}->{"qualdesc"}, $buffer2{$k}->{"qual"}, 2);
		} # End of for each statement
		close ($out2);
	} # End of if statemnet
} # End of if statement
else {
	# Print singletons
	my @keys = keys %buffer1;
	foreach my $k (@keys) {
		my $seqdesc = $buffer1{$k}->{"seqdesc"};
		my $qualdesc = $buffer1{$k}->{"qualdesc"};

		if ($buffer1{$k}->{"seqdesc"} =~ m/\/(1|2)$/) {
			$seqdesc = $`;
		} # end of if statement

		if ($buffer1{$k}->{"qualdesc"} =~ m/\/(1|2)$/) {
			$qualdesc = $`;
		} # end of if statement

		&printFastqSequence($out1, $seqdesc, $buffer1{$k}->{"seq"},
		                           $qualdesc, $buffer1{$k}->{"qual"}, 1);
		&printFastqSequence($out2, $seqdesc, "",
		                           $qualdesc, "", 2);
	} # end of for each statement
	
	@keys = keys %buffer2;
	foreach my $k (@keys) {
		my $seqdesc = $buffer2{$k}->{"seqdesc"};
		my $qualdesc = $buffer2{$k}->{"qualdesc"};
		
		if ($buffer2{$k}->{"seqdesc"} =~ m/\/(1|2)$/) {
			$seqdesc = $`;
		} # end of if statement

		if ($buffer2{$k}->{"qualdesc"} =~ m/\/(1|2)$/) {
			$qualdesc = $`;
		} # end of if statement
		
		&printFastqSequence($out1, $seqdesc, "",
		                           $qualdesc, "", 1);
		&printFastqSequence($out2, $seqdesc, $buffer2{$k}->{"seq"},
		                           $qualdesc, $buffer2{$k}->{"qual"}, 2);
	} # end of for each statement
	
	close ($out1);
	close ($out2);
} # end of else statement

