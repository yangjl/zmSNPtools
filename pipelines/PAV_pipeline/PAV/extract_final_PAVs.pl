#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

use constant true => 1;
use constant false => 0;

my ($fasta, @gff3, $output);
my $result = &GetOptions("fasta|f=s{1}" => \$fasta,
                         "gff3|g=s{1,}" => \@gff3,
						 "output|o=s{1}" => \$output);

unless ($result && defined($fasta) && scalar(@gff3) > 0 && defined($output)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl --fasta <fasta> --gff3 <gff3 files> --output <out>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # End of unless statemetn

print STDERR sprintf("\n");

# Reading gff3 files
my %discard;
foreach my $f (@gff3) {
	my $fh = new FileHandle();
	open ($fh, $f) or die("Cannot open file\n");
	print STDERR sprintf(" o Reading '%s' ... ", $f);
	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my @fields = split(/\t/, $_);
			my $id;
			if ($fields[$#fields] =~ m/Name=(\S+);Note/) {
				$id = $1;
			} # end of if statemetn
			else {
				die("Cannot determine ID\n");
			} # end of else statement

			$discard{$id} = true;
		} # end of if statement
	} # end of while loop
	close ($fh);
	print STDERR sprintf("DONE\n");
} # end of foreach statmenet

my $fh = new FileHandle();
open ($fh, $fasta) or die("Cannot open fasta file\n");

my $ofh = new FileHandle();
open ($ofh, sprintf(">%s", $output)) or die("Cannot create output file\n");
print STDERR sprintf(" o Dumping output ... ");
my ($name, $seq);
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($name) && !exists $discard{$name}) {
				&formatSequence($ofh, $name, $seq);
			} # end of if statemetn

			$name = $1;
			$seq = "";
		} # end of if statemetn
		else {
			$seq .= $_;
		} # end of else statemetn
	} # end of if statmeent
} # End of while loop
close ($fh);

if (defined($name) && !exists $discard{$name}) {
	&formatSequence($ofh, $name, $seq);
} # end of if statemetn
close ($ofh);
print STDERR sprintf("DONE\n");
print STDERR sprintf("\n");
