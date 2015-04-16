#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Schnablelab::Tools;

use constant QUALITIES_PER_LINE => 25;
use constant true => 1;
use constant false => 0;

sub cleanupScores {
	my $scores = $_[0];
	
	# Trim leading and trailing spaces
	$scores =~ s/^\s+//g;
	$scores =~ s/\s+$//g;

	# keep looping until each score value is separated by only 1 space
    while (index($scores, "  ") >= 0) {
        $scores =~ s/  / /g;
    } # End of while loop
        
    return $scores;
} # End of sub cleanupScores

sub formatQualities {
    my ($ofh, $name, @scores) = @_;

	if ($name !~ m/^>/) {
	    print $ofh sprintf(">%s\n", $name);
	} # End of if statement
	else {
		print $ofh sprintf("%s\n", $name);
	} # end of else statement

    while (scalar(@scores) > 0) {
        my @sub = splice(@scores, 0, QUALITIES_PER_LINE);
        @sub = map { sprintf("%2s", $_) } @sub;
        print $ofh sprintf("%s\n", join(" ", @sub));
    } # End of while loop
} # End of sub formatQualities


if (scalar(@ARGV) != 2) {
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("perl $0 <parsed info file> <lucy qualities file>\n");
	exit();
} # End of else statement

my %trimming;
my $fh = new FileHandle();
open ($fh, $ARGV[0]) or die("Cannot open parsed info file\n");
while (<$fh>) {
	chomp;
	if ($_ !~ m/^#/ && length($_) != 0) {
		my ($name, $length, $left, $right, $trimlength) = split(/\t/, $_);
		if (!exists $trimming{$name}) {
			$trimming{$name} = sprintf("%d,%d,%d", $left, $right, $trimlength);
		} # End of if statement
		else {
			die("*** Duplicate ID: $name\n");
		} # End of else statement
	} # End of if statement
} # End of while loop
close ($fh);

my ($name, $scores);
open ($fh, $ARGV[1]) or die("Cannot open lucy quality scores file\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/) {
			if (defined($name)) {
				if (exists $trimming{$name}) {
					my @qualities = split(/\s/, &cleanupScores($scores));
					my ($left, $right, $trimlength) = split(/,/, $trimming{$name});
					my @good_qualities = splice(@qualities, $left - 1, $trimlength);		# Extract middle elements in the array

					if (scalar(@good_qualities) == $trimlength) {
						&formatQualities(\*STDOUT, $name, @good_qualities);
					} # End of if statement
					else {
						my $scoreslength = scalar(@qualities);
						die("*** There was an error removing left/right scores\n$name $trimlength != $scoreslength\n");
					} # End of else statement
				} # end of if statmenet
			} # End of if statement

			$name = $1;
			undef($scores);
		} # end of if statement
		else {
			if (!defined($scores)) {
				$scores = $_;
			} # end of if statement
			else {
				$scores .= sprintf(" %s", $_);
			} # End of else statement
		} # End of else statement
	} # End of if statement
} # End of while loop
close ($fh);

if (defined($name)) {
	if (exists $trimming{$name}) {
		my @qualities = split(/\s/, &cleanupScores($scores));
		my ($left, $right, $trimlength) = split(/,/, $trimming{$name});
		my @good_qualities = splice(@qualities, $left - 1, $trimlength);		# Extract middle elements in the array

		if (scalar(@good_qualities) == $trimlength) {
			&formatQualities(\*STDOUT, $name, @good_qualities);			
		} # End of if statement
		else {
			my $scoreslength = scalar(@qualities);
			die("*** There was an error removing left/right scores\n$name $trimlength != $scoreslength\n");
		} # End of else statement
	} # end of if statmenet
} # End of if statement
