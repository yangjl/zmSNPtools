#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use Schnablelab::Tools;

use constant true => 1;
use constant false => 0;
use constant DEFAULT_MIN_IDENTITY => 0.95;
use constant DEFAULT_MIN_COVERAGE => 0.90;
use constant DEFAULT_MIN_ALIGN_LENGTH => 30;
use constant DEFAULT_PROCESSORS_COUNT => 24;

sub extract_tails {
	my ($inputFile, $outputFile) = @_;
	my $count = 0;

	my $fh = new FileHandle();
	open ($fh, $inputFile) or die("Cannot open fasta file\n");

	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s", $outputFile)) or die("Cannot create temporary file\n");

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
							&formatSequence($ofh, $id, $tmp[$i]);
							$count++;
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
				&formatSequence($ofh, $id, $tmp[$i]);
				$count++;
			} # End of if statement
			
			$index += length($tmp[$i]);	# Update index position
		} # end of for each statement
	} # end of if statement
	close ($ofh);

	my $command = sprintf("sequence_length.pl -f %s > %s.length", $outputFile, $outputFile);
	system($command);

	return $count;
} # End of sub extract_tails

sub prepareForNextSource {
	my ($fasta, $gff3) = @_;
	my %discard;
	my $aligned = 0;

	my $fh = new FileHandle();
	open ($fh, $gff3) or die("Cannot open gff3 file\n");
	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my @fields = split(/\t/, $_);
			my $id;
			if ($fields[$#fields] =~ m/Name=(\S+);Note/) {
				$id = $1;
				if (!exists $discard{$id}) {
					$discard{$id} = true;
					$aligned++;
				} # end of if statement
			} # end of if statement
		} # end of if statmenet
	} # End of while loop
	close ($fh);

	$fh = new FileHandle();
	open ($fh, $fasta) or die("Cannot open file\n");
	my ($name, $seq);

	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s.tmp", $fasta)) or die("Cannot create output file\n");

	while (<$fh>) {
		chomp;
		if (length($_) != 0) {
			if ($_ =~ m/^>(\S+)/) {
				if (defined($name) && !exists $discard{$name}) {
					&formatSequence($ofh, $name, $seq);
				} # end of if statmeent

				$name = $1;
				$seq = "";
			} # end of if statement
			else {
				$seq .= $_;
			} # end of else statement
		} # End of if statement
	} # End of while loop
	close ($fh);

	if (defined($name) && !exists $discard{$name}) {
		&formatSequence($ofh, $name, $seq);
	} # end of if statmeent
	close ($ofh);

	# Replace old file with new file for next round
	my $command = sprintf("mv %s.tmp %s", $fasta, $fasta);
	system($command);

	# Update length file
	$command = sprintf("sequence_length.pl -f %s > %s.length", $fasta, $fasta);
	system($command);

	return $aligned;
} # end of sub prepareForNextSource

sub generateConsensus {
	my ($inputFasta, $outputFasta, @alignments) = @_;
	my %contigs;

	foreach my $f (@alignments) {
		my $fh = new FileHandle();
		open ($fh, $f) or die("Cannot open file\n");
		while (<$fh>) {
			chomp;
			if (length($_) != 0 && $_ !~ m/^#/) {
				my @fields = split(/\t/, $_);

				if ($fields[$#fields] =~ m/Target=(\S+)_(\d+)-(\d+) (\d+) (\d+)$/) {
					my $id = $1;
					my $id_start = $2;
					my $id_end = $3;
					my $align_start = &min($4, $5);
					my $align_end = &max($4, $5);

					# Fix coordinates
					$align_start = $align_start - 1 + $id_start;
					$align_end = $align_end - 1 + $id_start;

					if (exists $contigs{$id}) {
						$contigs{$id} .= sprintf(";%s-%s", $align_start, $align_end);
					} # end of if statement
					else {
						$contigs{$id} = sprintf("%s-%s", $align_start, $align_end);
					} # End of else statement
				} # end of if statement
			} # end of if statement
		} # end of while loop
		close ($fh);
	} # end of for each statement

	my $fh = new FileHandle();
	open ($fh, $inputFasta) or die("Cannot open input fasta\n");

	my $ofh = new FileHandle();
	open ($ofh, sprintf(">%s", $outputFasta)) or die("Cannot create outptu file\n");
	my ($name, $seq);
	while (<$fh>) {
		chomp;
		if (length($_) != 0) {
			if ($_ =~ m/^>(\S+)/) {
				if (defined($name)) {
					if (exists $contigs{$name}) {
						foreach my $t (split(/;/, $contigs{$name})) {
							my ($s, $e) = split(/-/, $t);
							my $l = $e - $s + 1;
							my $m = "X" x $l;	# Create mask
							substr($seq, $s - 1, $l, $m);
						} # end of for each statmenet
						&formatSequence($ofh, $name, $seq);
					} # end of if statement
					else {	# No alignments, dump the whole thing
						&formatSequence($ofh, $name, $seq);
					} # end of else statement
				} # End of if statment

				$name = $1;
				$seq = "";
			} # end of if statement
			else {
				$seq .= $_;
			} # end of else statement
		} # end of if statement
	} # End of while loop
	close($fh);
				
	if (defined($name)) {
		if (exists $contigs{$name}) {
			foreach my $t (split(/;/, $contigs{$name})) {
				my ($s, $e) = split(/-/, $t);
				my $l = $e - $s + 1;
				my $m = "X" x $l;	# Create mask
				substr($seq, $s - 1, $l, $m);
			} # end of for each statmenet
			&formatSequence($ofh, $name, $seq);
		} # end of if statement
		else {	# No alignments, dump the whole thing
			&formatSequence($ofh, $name, $seq);
		} # end of else statement
	} # End of if statment
	close ($ofh);
} # End of generateConsensus

my ($fasta, $min_ident, $min_cvg, $min_align_length, $proc, $prefix, $dbpath);

my $result = &GetOptions("fasta|f=s{1}" => \$fasta,
                         "dbpath|db=s{1}" => \$dbpath,
                         "ident|i:f{1}" => \$min_ident,
						 "coverage|c:f{1}" => \$min_cvg,
						 "align|a:i{1}" => \$min_align_length,
						 "out|o=s{1}" => \$prefix,
						 "proc|p:i{1}" => \$proc);

unless ($result && defined($fasta) && defined($prefix) && defined($dbpath)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --fasta <fasta> --db <dbpath> --out <prefix> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS\n");
	print STDERR sprintf("  --ident <float>      : Specify minimum identity percentage required [DEFAULT: %2.2f]\n", DEFAULT_MIN_IDENTITY);
	print STDERR sprintf("  --coverage <float>   : Specify minimum coverage percentage required [DEFAULT: %2.2f]\n", DEFAULT_MIN_COVERAGE);
	print STDERR sprintf("  --align <int>        : Specify minimum number of bases required in the alignment [DEFAULT: %s]\n", DEFAULT_MIN_ALIGN_LENGTH);
	print STDERR sprintf("  --proc <int>         : Specify number of processors to use in alignment [DEFAULT: %s]\n", DEFAULT_PROCESSORS_COUNT);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

# Assign default values
my $run_id = $prefix;
$min_ident = DEFAULT_MIN_IDENTITY if (!defined($min_ident) || $min_ident !~ m/^\d+(\.\d+)?/);
$min_cvg = DEFAULT_MIN_COVERAGE if (!defined($min_cvg) || $min_cvg !~ m/^\d+(\.\d+)?/);
$min_align_length = DEFAULT_MIN_ALIGN_LENGTH if (!defined($min_align_length) || $min_align_length !~ m/^\d+$/);
$proc = DEFAULT_PROCESSORS_COUNT if (!defined($proc) || $proc !~ m/^\d+$/);

my $logFH = new FileHandle();
open ($logFH, sprintf(">%s.iterations.final.log", $prefix)) or die("Cannot create log file\n");

print $logFH sprintf("# %s\n", scalar(localtime(time)));
print $logFH sprintf("#\n");
print $logFH sprintf("# INPUT FASTA: %s\n", $fasta);
print $logFH sprintf("# MIN. IDENTITY: %s\n", $min_ident);
print $logFH sprintf("# MIN. COVERAGE: %s\n", $min_cvg);
print $logFH sprintf("# MIN. ALIGNMENT LENGTH: %s\n", $min_align_length);
print $logFH sprintf("# CPU COUNT: %s\n", $proc);

# Perform alignments
my $iter = 0;
my $total_aligned = 0;
do {
	$iter++;	# Increment
	$total_aligned = 0;		# Reset
	my $inputCount = 0;

	print $logFH sprintf("\n");
	print STDERR sprintf("\n");
	# Extract partial sequences for alignment
	print $logFH sprintf("--- ITERATION %s ---\n", &formatNumber($iter));
	print STDERR sprintf("--- ITERATION %s ---\n", &formatNumber($iter));
	print $logFH sprintf("  o Extracting partial sequences for alignment ... ");
	print STDERR sprintf("  o Extracting partial sequences for alignment ... ");
	
	my $count = &extract_tails($iter == 1 ? $fasta : sprintf("%s.iter%s.masked.fas", $run_id, $iter - 1), 
	                           sprintf("%s.iter%s.tails.fas", $run_id, $iter));
	$inputCount = $count;
	print $logFH sprintf("DONE [ %s sequences ]\n", &formatNumber($count));
	print STDERR sprintf("DONE [ %s sequences ]\n", &formatNumber($count));
	
	print $logFH sprintf("  o Performing tail re-alignemnts\n");
	print STDERR sprintf("  o Performing tail re-alignemnts\n");

	my @aligned_gff3;

	# AGPv3
	print $logFH sprintf("      + AGPv3 ... "); 
	print STDERR sprintf("      + AGPv3 ... "); 
	my $source = "AGPv3";
	my $db = "AGPv3+mito+chlo";
	my $command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                      sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	my $aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	my $remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration

	# AGPv2
	print $logFH sprintf("      + AGPv2 ... "); 
	print STDERR sprintf("      + AGPv2 ... "); 
	$source = "AGPv2";
	$db = "AGPv2+mito+chlo";
	$command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                   sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	$aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration

	# AGPv1
	print $logFH sprintf("      + AGPv1 ... "); 
	print STDERR sprintf("      + AGPv1 ... "); 
	$source = "AGPv1";
	$db="AGPv1+mito+chlo";
	$command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                   sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	$aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration
	
	# WUGSC
	print $logFH sprintf("      + WUGSC BACs ... "); 
	print STDERR sprintf("      + WUGSC BACs ... "); 
	$source = "WUGSC";
	$db = "WUGSC_BACs-20130313.pieces";
	$command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                   sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	$aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration
	
	# Nature Genetics
	print $logFH sprintf("      + Nature Genetics ... "); 
	print STDERR sprintf("      + Nature Genetics ... "); 
	$source = "Nature_Genetics";
	$db = "Nat_Genet_B73.novel";
	$command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                   sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	$aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration
	
	# MAGIv3.1
	print $logFH sprintf("      + MAGIv3.1 Contigs ... "); 
	print STDERR sprintf("      + MAGIv3.1 Contigs ... "); 
	$source = "MAGIv3.1";
	$db = "MAGIv3.1.contigs";
	$command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                   sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	$aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration
	
	# MAGIv4.0
	print $logFH sprintf("      + MAGIv4.0 Contigs ... "); 
	print STDERR sprintf("      + MAGIv4.0 Contigs ... "); 
	$source = "MAGIv4.0";
	$db = "MAGIv4.0.contigs";
	$command = sprintf("blastn -db %s -query %s.iter%s.tails.fas -out %s.iter%s.tails.%s.blastn.tabular -dust yes -num_threads %s -max_target_seqs 3 -outfmt 6 -task \"blastn\" > /dev/null 2>&1", 
	                   sprintf("%s/%s", $dbpath, $db), $run_id, $iter, $run_id, $iter, $source, $proc);
	system($command);
	$command = sprintf("blastn_tabular2gff3.pl -ql %s.iter%s.tails.fas.length -sl %s -blast %s.iter%s.tails.%s.blastn.tabular -o %s.iter%s.tails.%s.blastn.gff3 -it -i %s -c %s -a %s", 
	                   $run_id, $iter, sprintf("%s/%s.length", $dbpath, $db), $run_id, $iter, $source, $run_id, $iter, $source, $min_ident, $min_cvg, $min_align_length);
	system($command);
	unlink(sprintf("%s.iter%s.tails.%s.blastn.tabular", $run_id, $iter, $source));
	$aligned = &prepareForNextSource(sprintf("%s.iter%s.tails.fas", $run_id, $iter), sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$total_aligned += $aligned;
	push(@aligned_gff3, sprintf("%s.iter%s.tails.%s.blastn.gff3", $run_id, $iter, $source));
	$remaining = $count - $aligned;
	print $logFH sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	print STDERR sprintf("DONE [ %s / %s = %2.1f%% aligned ; %s / %s = %2.1f%% remaining ]\n",
	                     &formatNumber($aligned), &formatNumber($count), $count != 0 ? ($aligned / $count) * 100 : 0,
						 &formatNumber($remaining), &formatNumber($count), $count != 0 ? ($remaining / $count) * 100 : 0);
	$count = $remaining;	# Update count for next iteration

	# At this point we finished aligning all the tails to the B73 sources,
	# Generate a masked consensus for this iteration
	if ($iter == 1) {	# First iteration
		print $logFH sprintf("  o Masking '%s' ... ", $fasta);
		print STDERR sprintf("  o Masking '%s' ... ", $fasta);
		&generateConsensus($fasta, sprintf("%s.iter%s.masked.fas", $run_id, $iter), @aligned_gff3);
		print $logFH sprintf("DONE\n");
		print STDERR sprintf("DONE\n");
		print $logFH sprintf("  o Generated '%s' for next iteration\n", sprintf("%s.iter%s.masked.fas", $run_id, $iter));
		print STDERR sprintf("  o Generated '%s' for next iteration\n", sprintf("%s.iter%s.masked.fas", $run_id, $iter));
	} # End of if statement
	else {
		print $logFH sprintf("  o Masking '%s' ... ", sprintf("%s.iter%s.masked.fas", $run_id, $iter - 1));
		print STDERR sprintf("  o Masking '%s' ... ", sprintf("%s.iter%s.masked.fas", $run_id, $iter - 1));
		&generateConsensus(sprintf("%s.iter%s.masked.fas", $run_id, $iter - 1), sprintf("%s.iter%s.masked.fas", $run_id, $iter), @aligned_gff3);
		print $logFH sprintf("DONE\n");
		print STDERR sprintf("DONE\n");
		print $logFH sprintf("  o Generated '%s' for next iteration\n", sprintf("%s.iter%s.masked.fas", $run_id, $iter));
		print STDERR sprintf("  o Generated '%s' for next iteration\n", sprintf("%s.iter%s.masked.fas", $run_id, $iter));
	} # end of else statement

	unlink(sprintf("%s.iter%s.tails.fas", $run_id, $iter));
	unlink(sprintf("%s.iter%s.tails.fas.length", $run_id, $iter));

	# Clean up
	foreach my $f (@aligned_gff3) {
		unlink($f);
	} # End of for each statement

	print $logFH sprintf("  o Total partial tails aligned: %s / %s = %2.1f%%\n", &formatNumber($total_aligned), &formatNumber($inputCount),
	                     $inputCount != 0 ? ($total_aligned / $inputCount) * 100 : 0);
	print STDERR sprintf("  o Total partial tails aligned: %s / %s = %2.1f%%\n", &formatNumber($total_aligned), &formatNumber($inputCount),
	                     $inputCount != 0 ? ($total_aligned / $inputCount) * 100 : 0);
} while ($total_aligned != 0);

print $logFH sprintf("\n");
print STDERR sprintf("\n");
print $logFH sprintf("--- FINALIZING ---\n");
print STDERR sprintf("--- FINALIZING ---\n");
print $logFH sprintf("  o TOTAL ITERATIONS PERFORMED: %s\n", &formatNumber($iter - 1));
print STDERR sprintf("  o TOTAL ITERATIONS PERFORMED: %s\n", &formatNumber($iter - 1));
print $logFH sprintf("  o Renaming %s.iter%s.masked.fas as final output %s.itererations.final.fas ... ",
                     $run_id, $iter, $run_id);
print STDERR sprintf("  o Renaming %s.iter%s.masked.fas as final output %s.itererations.final.fas ... ",
                     $run_id, $iter, $run_id);
my $command = sprintf("mv %s.iter%s.masked.fas %s.iterations.final.fas", $run_id, $iter, $run_id);
system($command);
print $logFH sprintf("DONE\n");
print STDERR sprintf("DONE\n");

print $logFH sprintf("  o Cleanup ... ");
print STDERR sprintf("  o Cleanup ... ");
for (my $i=0; $i < $iter; $i++) {
	unlink(sprintf("%s.iter%s.masked.fas", $run_id, $i));
} # end of for loop
print $logFH sprintf("DONE\n");
print STDERR sprintf("DONE\n");

print $logFH sprintf("\n");
print STDERR sprintf("\n");

close ($logFH);
