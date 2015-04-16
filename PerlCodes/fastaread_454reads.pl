#!/usr/bin/perl -w
# fastaread_454reads.pl 
# Usage: perl fastaread_454reads.pl --i <sequence name list> --d <454 fasta file or trimmed sequence file> 

use strict;
use warnings;
use Getopt::Long;

my $name;
my $seq;
my %name_seq;
my ($input, $db, $help);

GetOptions ("i=s" => \$input, "d=s" => \$db, "help" => \$help);
if ($help) {
	print <<EOF
	Usage: perl fastaread_454reads.pl --i <sequence name list> --d <454 fasta file or trimmed sequence file>
	the script automatically search the fixed format of 454 sequence name in the input file and look up the 
	corresponding sequence in the 454 database;
EOF
}

# open the db file:
open(IN, $db) || die "Cannot open db file: $db\n";
while (<IN>) {
	chomp;
	if (/^>(\d+_\d+_\d+)/) {
		if (defined $name) {
			$name_seq{$name} = $seq;
		}
		$name = $1;
		$seq = '';
	} else {
		$seq .= $_;
	}
}
$name_seq{$name} = $seq;
close IN;

# open the input file:
open(IN, $input) || die "Cannot open $input!";
while (<IN>) {
	$_ =~ /(\d+_\d+_\d+)/;
	if (exists $name_seq{$1}) {
		print ">$1\n$name_seq{$1}\n";
	} else {
		print STDERR "Cannot find sequence: $1.\n";
	}
}
close IN;
