#!/usr/bin/perl -w

# blast.filter.pl
# usage: perl blast.filter.pl <blastnparser.pl output> [--identity int] [--coverage int] [--subCov int] [--len int] [--qStart int] [--qEnd int] [--subStart int] [--subEnd int] [--e int] [--bit int] [--tail int]
use strict;
use warnings;
use Getopt::Long;

sub errInf {
	print <<EOF;
	Usage: perl blast.filter.pl <blastnparser.pl output> [--identity int] [--coverage int] [--subCov int] [--len int] [--qStart int] [--qEnd int] [--subStart int] [--subEnd int] [--e int] [--bit int] [--tail int]
	This is the filter for blastn parsing file (output of blastnparser.pl).
	The output should have: 
		>= identity (identity <=100)
		>= coverage (100*overlap length/length of query, <=100)
		>= subCov (100*overlap length/length of subject, <=100)
		>= len (overlap length)
		<= qStart (query start)
		<= qEnd (query end)
		<= subStart (subject start)
		<= subEnd (subject end)
		<= e (e-value, the input value is the factorial number, e.g. 15 represent "e-15", default e = 0)
		>= bit (bit score) 
		<= tail (the total overhang at both sides)
		<= subPosition (the minimum number of Sub.Start and Sub.End is smaller than "subPosition")
EOF
}

# parameter variables:
my $identity;
my $len=0;
my $coverage;
my $subCov;
my $qStart; 
my $qEnd;
my $subStart;
my $subEnd;
my $e=0;
my $bit=0;
my $tail; # the total overhang at both sides
my $subPos; # the minimum number of Sub.Start and Sub.End is smaller thant "subPosition"

# other variables:
my @line;
my $file;

# read the parameter: 
my $result = &GetOptions("identity=i" => \$identity, "coverage=i" => \$coverage, "subCov=i" => \$subCov, "len=i" => \$len, "qStart=i" => \$qStart, "qEnd=i" => \$qEnd, "subStart=i" => \$subStart, "subEnd=i" => \$subEnd, "e=i" => \$e, "bit=i" => \$bit, "tail=i" => \$tail, "subPosition=i" => \$subPos);
if (!$result) {
	&errInf;
	exit;
}

open (IN, $ARGV[0]) || die "\tCannot open the input file: $ARGV[0].\n";# $file = 0;
# if (!$file) {
#	print "\tCannot open the input file: $ARGV[0].\n";
#	&errInf;
#	exit;
#}
my $evalue = 10**(-$e);
my $qStart_value;
my $qEnd_value;
my $subStart_value;
my $subEnd_value;
my $identity_value;
my $coverage_value;
my $subCov_value;
my $bit_value;
my $tail_value; 
my $subPos_value;

# skip the 1st line, headline:
$_ = <IN>;

# read the input file from the 2nd line
while (<IN>) {
	chomp;
	# if "No hits found", skip this line:
	unless (/No hits found/) {
		@line= split(/\t/,$_);
		$line[9] =~ /(\d+\.\d+)\%/;
		$coverage_value = $1;
		$subCov_value = abs($line[17]-$line[18])*100/$line[6]; 
		$line[11] =~ /(\d+)\%/;
		$identity_value  = $1;
		$line[7] =~ /(\d+)\sbits/;
		$bit_value = $1;
		# put the minumum value of subStart and subEnd into $subPos_value:
		my @sub=();
		@sub = sort {$a <=> $b} ($line[17],$line[18]); 
		$subPos_value = $sub[0];
	
		# compare values in each line to the input parameter
		if ((defined $identity) && ($identity_value < $identity)) {
			$identity_value = 0;
		} else {$identity_value = 1;}
	
		if ((defined $coverage) && ($coverage_value < $coverage)) {
			$coverage_value = 0;
		} else {$coverage_value = 1;}
		
		if ((defined $subCov) && ($subCov_value < $subCov)) {
			$subCov_value = 0;
		} else {$subCov_value = 1;}
		
		if ((defined $qStart) && ($line[15] > $qStart)) {
			$qStart_value = 0;
		} else {$qStart_value = 1;}
	
		if ((defined $qEnd) && ($line[16] > $qEnd)) {
			$qEnd_value = 0;
		} else {$qEnd_value = 1;}
	
		if ((defined $subStart) && ($line[17] > $subStart)) {
			$subStart_value = 0;
		} else {$subStart_value = 1;}
	
		if ((defined $subEnd) && ($line[18] > $subEnd)) {
			$subEnd_value = 0;
		} else {$subEnd_value = 1;}
	
		if ((defined $bit) && ($bit_value < $bit)) {
			$bit_value = 0;
		} else {$bit_value = 1;}
		
		# subPos vs $subPos_value ( the minimum of subStart and subEnd)
		if ((defined $subPos) && ($subPos_value > $subPos)) {
			$subPos_value = 0;
		} else {$subPos_value = 1;}
	
		if ((defined $tail) && (($line[2]-$line[16]+$line[15]-1) > $tail)) {
			$tail_value = 0;
		} else {$tail_value = 1;}
	
		# output all lines satisfied to the criteria:
		if ($identity_value && $coverage_value && (($line[16]-$line[15]) >= $len) && $subCov_value && $qStart_value && $qEnd_value && $subStart_value && $subEnd_value && ($line[8]<=$evalue) && $bit_value && $tail_value && $subPos_value) {
			print "$_\n";
		}
	}
}
close IN;
