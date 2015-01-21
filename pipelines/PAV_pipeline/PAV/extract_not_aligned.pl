#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

use constant true => 1;
use constant false => 0;

sub getName {
	my @tokens = split(/;/, $_[0]);
	my $found = false;
	my $name;

	for (my $i=0; $i < scalar(@tokens) && !$found; $i++) {
		if ($tokens[$i] =~ m/^Name=/) {
			$name = $';
			$found = true;
		} # end of if statement
	} # end of for loop

	return $name;
} # End of getName

my ($gff3, $fasta, $output);
my $result = &GetOptions("gff3|g=s{1}" => \$gff3,
                         "fasta|f=s{1}" => \$fasta,
						 "output|o=s{1}" => \$output);

unless ($result && defined($gff3) && defined($fasta) && defined($output)) {
	print STDERR sprintf("perl %s --gff3 <gff3 alignments> --fasta <fasta file> --output <out>\n", $0);
	exit();
} # End of unless statement

# Reading alignments and keep track of aligned names
my %aligned;
my $fh = new FileHandle();
open ($fh, $gff3) or die("Cannot open gff3 file\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0 && $_ !~ m/^#/) {
		my @fields = split(/\t/, $_);
		if ($fields[2] !~ m/^HSP$/i) {
			my $name = &getName($fields[$#fields]);
			$aligned{$name} = true;
		} # end of if statement
	} # end of if statement
} # End of while loop
close ($fh);

# Reading fasta file
$fh = new FileHandle();
open ($fh, $fasta) or die("Cannot open fasta file\n");
my ($name, $seq);

my $ofh = new FileHandle();
open ($ofh, sprintf(">%s", $output)) or die("Cannot create output fasta file\n");
my $count = 0;
my $new= 0;
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($name) && !exists $aligned{$name}) {
				&formatSequence($ofh, $name, $seq);
				$new++;
			} # End of if statement

			$name = $1;
			$seq = "";
			$count++;
		} # end of if statement
		else {
			$seq .= $_;
		} # end of else statement
	} # end of if statement
} # End of while loop

if (defined($name) && !exists $aligned{$name}) {
	&formatSequence($ofh, $name, $seq);
	$new++;
} # End of if statement

close ($fh);
close ($ofh);
