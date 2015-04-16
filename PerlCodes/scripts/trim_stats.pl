#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

sub removeComma {
	my $str = $_[0];
	$str =~ s/,//g;

	return $str;
} # End of sub removeComma

use constant LAST_N_LINES => 20;	# parse the last 20 lines for each file

my @log_files;
my $output;
my $result = &GetOptions("log|l|files|f=s{1,}" => \@log_files, "output|o=s{1}" => \$output);

unless ($result && scalar(@log_files) > 0 && defined($output)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --log|-l|--files|-f <*.log or *.log.gz files> --output|-o <output>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

print STDERR sprintf("\n");

my %values;
foreach my $f (@log_files) {
	print STDERR sprintf(" o Processing '%s' ... ", $f);
	my $logFH = new FileHandle();

	if ($f =~ m/\.gz$/) {
		open ($logFH, sprintf("zcat %s | tail -n %s|", $f, LAST_N_LINES)) or 
			die(sprintf("cannot parse last %s lines of %s\n", LAST_N_LINES, $f));
	} # end of if statemnt
	else {
		open ($logFH, sprintf("tail -n %s %s|", LAST_N_LINES, $f)) or 
			die(sprintf("cannot parse last %s lines of %s\n", LAST_N_LINES, $f));
	} # end of else statemetn

	my ($raw, $trimmed, $raw_bp, $trimmed_bp, $raw_avg, $trimmed_avg);
	while (<$logFH>) {
		chomp;
		if ($_ =~ m/^# TOTAL READS \(RAW\): (\S+)/) {
			$raw = &removeComma($1);
		} # End of if statement
		elsif ($_ =~ m/^# TOTAL READS REMAINING \(POST-TRIM\): (\S+)/) {
			$trimmed = &removeComma($1);
		} # end of else if statmenet
		elsif ($_ =~ m/^# TOTAL BASE PAIRS \(RAW\): (\S+)/) {
			$raw_bp = &removeComma($1);
		} # end of else if statement
		elsif ($_ =~ m/^# TOTAL BASE PAIRS REMAINING \(POST-TRIM\): (\S+)/) {
			$trimmed_bp = &removeComma($1);
		} # end of if statement
		elsif ($_ =~ m/^# AVERAGE READ LENGTH \(RAW\): (\S+)/) {
			$raw_avg = &removeComma($1);
		} # end of else if statement
		elsif ($_ =~ m/^# AVERAGE READ LENGTH \(POST-TRIM\): (\S+)/) {
			$trimmed_avg = &removeComma($1);
		} # end of else if statmenet
	} # end of while loop
	close ($logFH);

	if (defined($raw) && defined($trimmed) && defined($raw_bp) && defined($trimmed_bp) &&
	    defined($raw_avg) && defined($trimmed_avg)) {
		my @tmp = split(/\//, $f);
		my $name = pop(@tmp);
		my ($prefix, $mate);
		if ($name =~ m/-(1|2)\.trimmed\.log(\.gz)?$/ || $name =~ m/\.trimmed-(1|2)\.log(\.gz)?$/) {
			$prefix = $`;
			$mate = $1;
		} # end of if statement
		else {
			$prefix = $f;
			$mate = 1;
		} # end of else statment

		# Remove *.trimmed.log or *.trimmed.log.gz
		$name =~ s/\.trimmed.log(\.gz)?$//g;

		$values{$prefix}->{$mate}->{"filename"} = $name;
		$values{$prefix}->{$mate}->{"raw"} = $raw;
		$values{$prefix}->{$mate}->{"trimmed"} = $trimmed;
		$values{$prefix}->{$mate}->{"raw_bp"} = $raw_bp;
		$values{$prefix}->{$mate}->{"trimmed_bp"} = $trimmed_bp;
		$values{$prefix}->{$mate}->{"raw_avg"} = $raw_avg;
		$values{$prefix}->{$mate}->{"trimmed_avg"} = $trimmed_avg;
		print STDERR sprintf("DONE\n");
	} # end of if statment
	else {
		print STDERR sprintf("FAILED\n");
		print STDERR sprintf("\nThere were errors determining values in the specified log file\n");
		print STDERR sprintf("\n");
		exit();
	} # end of else statmenet
} # end of for each statement

# Dumping output
my @keys = sort {$a cmp $b} keys %values;
if (scalar(@keys) > 0) {
	print STDERR sprintf(" o Dumping output ... ");
	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s", $output)) or die("Cannot create output file\n");
	print $ofh sprintf("# File\t No. Raw Reads\tRaw Base Pairs (bp)\tRead Length (bp)\t");
	print $ofh sprintf("No. Trimmed Reads\t%% Trimmed Reads Remaining\tTrimmed Base Pairs (bp)\t%% Trimmed Bases Remaining\tRead Length (bp)\t");
	print $ofh sprintf("No. Trimmed Reads\tTrimmed Base Pairs (bp)\tRead Length (bp)\n");

	foreach my $k (@keys) {
		my @mates = sort {$a <=> $b} keys %{ $values{$k} };
		my ($total_raw, $total_trimmed, $total_raw_bp, $total_trimmed_bp, $total_raw_avg, $total_trimmed_avg) = (0, 0, 0, 0, 0, 0);

		foreach my $m (@mates) {
			my $name;
			if (scalar(@mates) == 1) {
				$name = $values{$k}->{$m}->{"filename"};
			} # end of if statemetn
			else {
				$name = sprintf("%s-%s", $k, $m);
			} # end of if statement

			my $raw = $values{$k}->{$m}->{"raw"}; $total_raw += $raw;
			my $trimmed = $values{$k}->{$m}->{"trimmed"}; $total_trimmed += $trimmed;
			my $raw_bp = $values{$k}->{$m}->{"raw_bp"}; $total_raw_bp += $raw_bp;
			my $trimmed_bp = $values{$k}->{$m}->{"trimmed_bp"}; $total_trimmed_bp += $trimmed_bp;
			my $raw_avg = $values{$k}->{$m}->{"raw_avg"}; $total_raw_avg += $raw_avg;
			my $trimmed_avg = $values{$k}->{$m}->{"trimmed_avg"}; $total_trimmed_avg += $trimmed_avg;

			my $trim_percentage = ($trimmed / $raw) * 100;
			my $trim_bp_percentage = ($trimmed_bp / $raw_bp) * 100;

			print $ofh sprintf("%s\t%s\t%s\t%s\t", $name, &formatNumber($raw), &formatNumber($raw_bp), &formatNumber($raw_avg));
			print $ofh sprintf("%s\t%2.1f%%\t%s\t%2.1f%%\t%s\t", &formatNumber($trimmed), $trim_percentage,
																 &formatNumber($trimmed_bp), $trim_bp_percentage,
																 &formatNumber($trimmed_avg));
			print $ofh sprintf("%s (%2.1f%%)\t%s (%2.1f%%)\t%s\n", &formatNumber($trimmed), $trim_percentage,
																 &formatNumber($trimmed_bp), $trim_bp_percentage,
																 &formatNumber($trimmed_avg));
		} # end of for each statement

		if (scalar(@mates) > 1) {
			my $total_trim_percentage = ($total_trimmed / $total_raw) * 100;
			my $total_trim_bp_percentage = ($total_trimmed_bp / $total_raw_bp) * 100;
			
			print $ofh sprintf("%s total\t%s\t%s\t%s\t", $k, &formatNumber($total_raw), &formatNumber($total_raw_bp), &formatNumber(int($total_raw_avg / scalar(@mates))));
			print $ofh sprintf("%s\t%2.1f%%\t%s\t%2.1f%%\t%s\t", &formatNumber($total_trimmed), $total_trim_percentage,
																 &formatNumber($total_trimmed_bp), $total_trim_bp_percentage,
																 &formatNumber(int($total_trimmed_avg / scalar(@mates))));
			print $ofh sprintf("%s (%2.1f%%)\t%s (%2.1f%%)\t%s\n", &formatNumber($total_trimmed), $total_trim_percentage,
																 &formatNumber($total_trimmed_bp), $total_trim_bp_percentage,
																 &formatNumber(int($total_trimmed_avg / scalar(@mates))));
			print $ofh sprintf("\n");
		} # End of if statement
	} # end of for each statement
	close ($ofh);

	print STDERR sprintf("DONE\n");
} # end of if statmenet

print STDERR sprintf("\n");
