#!/usr/bin/perl -w
# fastaread_454reads.pl 
# Usage: perl fastaread_454reads.pl --i <sequence name list> --d <454 fasta file or trimmed sequence file> 

use strict;
use warnings;
use Getopt::Long;

use constant col => 6;

my @line;
my $name;
my $seq;
my %name_seq;
my ($input, $db);

GetOptions ("i=s" => \$input, "d=s" => \$db);

# open the input file:
open(IN, $input) || die "Cannot open $input!";
while (<IN>) {
	if (/(\d+_\d+_\d+)/) {
		$454{$1}++;
	}
}
close IN;


# open the db file:
open(IN, $db) || die "Cannot open db file: $db\n";
while (<IN>) {
	chomp;
	if (/^>(\d+_\d+_\d+)/) {
		if (defined $name) {
			if (exists $454{$name}) {
				$name_seq{$name} = $seq;
			}
		}
		$name = $1;
		$seq = '';
	} else {
		$seq .= $_;
	}
}
if (exists $454{$name}) {
	$name_seq{$name} = $seq;
}
close IN;

my %found;
# open the input file:
open(IN, $input) || die "Cannot open $input!";
while (<IN>) {
	@line = split(/\t/,$_);
	$_ =~ /(\d+_\d+_\d+)/;
	if (exists $name_seq{$1}) {
		if (!exists $found{$line[1].$1}) {
			print ">$line[1]_$1\n$name_seq{$1}\n";
			$found{$line[1].$1}++;
		}
	} else {
		print STDERR "Cannot find sequence: $1.\n";
	}
}
close IN;

