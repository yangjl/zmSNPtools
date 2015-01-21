#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use constant DEFAULT_TYPE => "fasta";

my @files;
my $type;
my $result = &GetOptions("files|f=s{1,}" => \@files, "type|t:s{1}" => \$type);

unless ($result && scalar(@files) > 0) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --files|-f <files> [--type <fasta|fastq>]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("Note: Default \"type\" is %s\n", DEFAULT_TYPE);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

$type = DEFAULT_TYPE if (!defined($type));

if (lc($type) eq "fasta" || lc($type) eq "fastq") {
	if (lc($type) eq "fasta") {
		foreach my $f (@files) {
			my ($name, $seq);
			my $fh = new FileHandle();
			open ($fh, $f) or die("Cannot open fasta file\n");
			while (<$fh>) {
				chomp;
				if (length($_) != 0) {
					if ($_ =~ m/^>(\S+)/) {
						if (defined($name)) {
							printf("%s\t%s\n", $name, length($seq));
						} # End of if statement

						$name = $1;
						$seq = "";
					} # End of if statement
					else {
						$seq .= $_;
					} # end of else statemnet
				} # End of if statemnet
			} # end of while loop
			close ($fh);

			# very last sequence
			if (defined($name)) {
				printf("%s\t%s\n", $name, length($seq));
			} # End of if statement
		} # end of for each statement
	} # end of if statmenet
	else {
		foreach my $f (@files) {
			my $fh = new FileHandle();
			open ($fh, $f) or die("Cannot open fastq file\n");
			my ($seq_desc, $seq, $qual_desc, $qual);
			while (!eof($fh)) {
				$seq_desc = <$fh>;		chomp($seq_desc);
				$seq = <$fh>;			chomp($seq);
				$qual_desc = <$fh>;		chomp($qual_desc);
				$qual = <$fh>;			chomp($qual);

				$seq_desc = substr($seq_desc, 1);
				printf("%s\t%s\n", $seq_desc, length($seq));
			} # end of while loop
			close ($fh);
		} # end of foreach statement
	} # end of ele statement
} # end of if statement
else {
	print STDERR sprintf("\n");
	print STDERR sprintf("INVALID TYPE: %s\n", $type);
	print STDERR sprintf("\n");
} # end of else statement
