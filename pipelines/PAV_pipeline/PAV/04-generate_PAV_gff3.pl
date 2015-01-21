#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

my ($fasta, $source);
my $result = &GetOptions("fasta|f=s{1}" => \$fasta, "source|s=s{1}" => \$source);

unless ($result && defined($fasta) && defined($source)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --fasta <masked fasta> --source <str>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

printf("## gff3-version\t3\n");
printf("# %s\n", scalar(localtime(time)));
printf("#\n");
printf("# INPUT FASTA: %s\n", $fasta);
printf("# SOURCE: %s\n", $source);

my $fh = new FileHandle();
open ($fh, $fasta) or die("cannot open file\n");
my ($name, $seq);
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($name)) {
				if (&countNucleotides($seq, "X") != length($seq)) {
					printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s;Name=%s;Note=Length: %s bp\n",
						   $name, $source, "Contig", 1, length($seq), ".", ".", ".", $name, $name, length($seq));
					
					my $start = 1;
					my $end;
					my @tmp = split(/(X{1,})/, $seq);
					for (my $i=0; $i < scalar(@tmp); $i++) {
						if (length($tmp[$i]) != 0) {
							my $type;
							if ($tmp[$i] =~ m/^X{1,}/) {	# B73 Region
								$type = "B73-Like";
							} # End of if statement
							else {
								$type = "PAV";
							} # end of ele statmeent

							$end = $start + length($tmp[$i]) - 1;
							printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tParent=%s;ID=%s_%s-%s;Name=%s_%s-%s;Note=Length: %s bp\n",
								   $name, $source, $type, $start, $end, ".", ".", ".", $name, $name, $start, $end, $name, $start, $end,
								   $end - $start + 1);
							$start += length($tmp[$i]);
						} # end of if statement
					} # end of for loop
				} # end of if statemetn
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
	if (&countNucleotides($seq, "X") != length($seq)) {
		printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s;Name=%s;Note=Length: %s bp\n",
			   $name, $source, "Contig", 1, length($seq), ".", ".", ".", $name, $name, length($seq));
		
		my $start = 1;
		my $end;
		my @tmp = split(/(X{1,})/, $seq);
		for (my $i=0; $i < scalar(@tmp); $i++) {
			my $type;
			if ($tmp[$i] !~ m/^X{1,}/) {	# B73 Region
				$type = "B73-Like";
			} # End of if statement
			else {
				$type = "PAV";
			} # end of ele statmeent

			$end = $start + length($tmp[$i]) - 1;
			printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tName=%s_%s-%s;ID=%s_%s-%s;Parent=%s;Note=Length: %s bp;Target=%s %s %s\n",
				   $name, $source, $type, $start, $end, ".", ".", ".", $name, $start, $end, $name, $start, $end, $name,
				   $end - $start + 1, $name, $start, $end);
			$start += length($tmp[$i]);
		} # end of for loop
	} # end of if statement
} # end of if statement

