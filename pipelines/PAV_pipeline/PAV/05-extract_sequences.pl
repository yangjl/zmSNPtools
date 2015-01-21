#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use constant true => 1;
use constant false => 0;

my ($gff3, $file);
my $result = &GetOptions("gff3|g=s{1}" => \$gff3, "fasta|f=s{1}" => \$file);

unless ($result && defined($gff3) && defined($file)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --gff3 <gff3 file> --fasta <fasta file>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

my %keep;
my $fh = new FileHandle();
open ($fh, $gff3) or die("cannot open gff3 file\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0 && $_ !~ m/^#/) {
		my @fields = split(/\t/, $_);
		if ($fields[2] =~ m/^Contig$/ && $fields[$#fields] =~ m/Name=(\S+);Note/) {
			$keep{$1} = true;
		} # End of if statemnet
	} # End of if statemnet
} # End of while loop
close ($fh);

$fh = new FileHandle();
open ($fh, $file) or die("Cannot open fasta or quality file\n");
my $print = false;
while (<$fh>) {
	chomp;
	if ($_ =~ m/^>(\S+)/) {
		if (exists $keep{$1}) {
			$print = true;
		} # end of if statement
		else {
			$print = false;
		} # End of else statement
	} # End of if statement

	if ($print) {
		printf("%s\n", $_);
	} # End of if statement
} # End of while loop
close ($fh);

