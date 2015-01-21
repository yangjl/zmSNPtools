#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use Schnablelab::Tools;

my ($fasta);
my $result = &GetOptions("fasta|f=s{1}" => \$fasta);

unless ($result && defined($fasta)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --fasta <B73 masked fasta>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # end of unles statement

my $fh = new FileHandle();
open ($fh, $fasta) or die("Cannot open file\n");

my ($name, $seq);
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($name)) {
				my $index = 1;
				my @tmp = split(/(X{1,})/, $seq);
				for (my $i=0; $i < scalar(@tmp); $i++) {
					if (length($tmp[$i]) != 0 && $tmp[$i] !~ m/^X{1,}$/) {
						my $id = sprintf("%s_%s-%s", $name, $index, $index + length($tmp[$i]) - 1);
						&formatSequence(\*STDOUT, $id, $tmp[$i]);
					} # End of if statement
					
					$index += length($tmp[$i]);	# Update index position
				} # end of for each statement
			} # end of if statement

			$name = $1;
			$seq = "";
		} # end of if statement
		else {
			$seq .= $_;
		} # end of else statement
	} # end of if statemnet
} # end of while loop
close ($fh);
			
if (defined($name)) {
	my $index = 1;
	my @tmp = split(/(X{1,})/, $seq);
	for (my $i=0; $i < scalar(@tmp); $i++) {
		if ($tmp[$i] !~ m/^X{1,}$/) {
			my $id = sprintf("%s_%s-%s", $name, $index, $index + length($tmp[$i]) - 1);
			&formatSequence(\*STDOUT, $id, $tmp[$i]);
		} # End of if statement
		
		$index += length($tmp[$i]);	# Update index position
	} # end of for each statement
} # end of if statement
close ($fh);

