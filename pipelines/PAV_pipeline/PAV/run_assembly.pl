#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Time::Local;

use Schnablelab::Tools;

use constant true => 1;
use constant false => 0;

# Default user provided parameters
use constant DEFAULT_K => 85;
use constant DEFAULT_PROC => 24;

# Default ABySS parameters (DO NOT CHANGE)
use constant DEFAULT_ALIGNER => "map";
use constant DEFAULT_C => 10;
use constant DEFAULT_S => 400;
use constant DEFAULT_N => 10;
use constant DEFAULT_M => 80;
use constant DEFAULT_P => 0.95;

# Format seconds notation into hours:mins:seconds
sub formatTime {
    my $totaltime = $_[0];
    my $str;

    $str = sprintf("%02d:", $totaltime / 3600);    # total hours
    $totaltime = $totaltime % 3600;
    $str .= sprintf("%02d:", $totaltime / 60);      # total minutes
    $totaltime = $totaltime % 60;
    $str .= sprintf("%02d", $totaltime);            # total sconds
    return $str;
} # End of sub formatTime

# Detect processors
sub detectProcessorsCount {
	open (COUNT, sprintf("cat /proc/cpuinfo | grep -c \"processor\"|"));
	my $count = <COUNT>;	chomp($count);
	if ($count !~ m/^\d+$/) {
		$count = DEFAULT_PROC;
	} # end of if statement
	
	return $count;
} # end of sub detectProcessorsCount

my ($list, $inputDir, $outputDir, $k, $proc);

my $result = &GetOptions("list|l=s{1}" => \$list,
                         "inputDir|i=s{1}" => \$inputDir,
						 "outdir|o=s{1}" => \$outputDir,
						 "kmer|k:i{1}" => \$k,
						 "proc|p:i{1}" => \$proc);

# Set default processors
$proc = &detectProcessorsCount();

unless ($result && defined($list) && defined($inputDir) && defined($outputDir)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("     perl %s --list <genotyping_info.txt> --inputdir <indir> --outdir <outdir> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --list <genotyping_info.txt>          : File path containing genotyping info\n");
	print STDERR sprintf("   --inputdir <input dir>                : Directory containing gunzip repeatmasked fastq files for assembly\n");
	print STDERR sprintf("   --outdir <output dir>                 : Directory in which assembly results will be saved\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("   --kmer|-k <int>                       : K-Mer length to run assembly [ DFDAULT: %s ]\n", DEFAULT_K);
	print STDERR sprintf("   --proc|-p <int>                       : Number of processors to use for assembly [ DEFAULT: %s ]\n", $proc);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

# Assign default parameters
$k = DEFAULT_K if (!defined($k) || $k !~ m/^\d+$/);

# Remove extra '/' from input and output dirs
$inputDir = substr($inputDir, 0, length($inputDir) - 1) if ($inputDir =~ m/\//);
$outputDir = substr($outputDir, 0, length($outputDir) - 1) if ($outputDir =~ m/\//);

print STDERR sprintf("\n");

# Reading genotyping file
my %genotypes;
print STDERR sprintf(" o Reading genotype info file '%s' ... ", $list);
my $fh = new FileHandle();
open ($fh, $list) or die("Cannot open genotype info file\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0 && $_ !~ m/^#/) {
		my ($lib, $geno) = split(/\t/, $_);
		$genotypes{$geno}->{"lib"} = $lib;		# Saving info
	} # end of if statement
} # end of while loop
close ($fh);

my @genotypes = sort { $a cmp $b } keys %genotypes;
print STDERR sprintf("DONE (%s found)\n", &formatNumber(scalar(@genotypes)));

# Obtain list of compressed files
print STDERR sprintf(" o Obtaining list of compressed files for processing in '%s' ... ", $inputDir);
my $dirFH = new FileHandle();
opendir ($dirFH, $inputDir) or die("Cannot read input directory files\n");
my $filesCount = 0;
foreach my $f (readdir($dirFH)) {
	if ($f =~ m/^(\S+)\.(rep\d+)\.(\S+)\.trimmed\.repmasked-(paired|singletons)-(1|2)\.fq.gz$/) {
		my $geno = $1;
		my $rep = $2;
		my $lib = $3;
		my $type = $4;
		my $num = $5;

		# Saving file info
		if (exists $genotypes{$geno} && $genotypes{$geno}->{"lib"} eq $lib) {	# Matched file name
			if (!exists $genotypes{$geno}->{"reps"}->{$rep}->{$type}->{$num}) {
				$genotypes{$geno}->{"reps"}->{$rep}->{$type}->{$num} = $f;
				$filesCount++;
			} # end of if statement
		} # end of if statement
	} # end of if statement
} # end of for each statement
close ($dirFH);
print STDERR sprintf("DONE (%s files found)\n", &formatNumber($filesCount));

# Dummy variable, for formatting purposes
my $dummy = length(&formatNumber(scalar(@genotypes)));

# Checking files, generating assembly scripts, and running
for (my $ii=0; $ii < scalar(@genotypes); $ii++) {
	my $geno = $genotypes[$ii];
	print STDERR sprintf(" o [ %\Q$dummy\Es / %\Q$dummy\Es ]: '%s' Assembly\n", &formatNumber($ii+1), &formatNumber(scalar(@genotypes)), $geno);
	my $script_name = sprintf("%s.%s.k%s-%s.sh", $geno, $genotypes{$geno}->{"lib"}, $k, DEFAULT_ALIGNER);
	
	if (-d sprintf("%s/%s", $outputDir, $geno)) {
		print STDERR sprintf("    + IGNORED: %s/%s already present\n", $outputDir, $geno);
	} # end of if statemnet
	else {
		my @reps = sort { $a cmp $b } keys %{ $genotypes{$geno}->{"reps"} };
		print STDERR sprintf("    + Start Time: %s\n", scalar(localtime(time)));
		my $start_time = timelocal(localtime(time));
		my $ok = true;
		foreach my $r (@reps) {
			print STDERR sprintf("    + Checking '%s' files\n", $r);
			print STDERR sprintf("        - %s.%s.%s.trimmed.repmasked-paired-1.fq.gz ... ", $geno, $r, $genotypes{$geno}->{"lib"});
			if (!exists $genotypes{$geno}->{"reps"}->{$r}->{"paired"}->{1}) {
				print STDERR sprintf("NOT PRESENT\n");
				$ok = false;
			} # end of if statement
			else {
				print STDERR sprintf("PRESENT\n");
			} # End of else statement
			
			print STDERR sprintf("        - %s.%s.%s.trimmed.repmasked-paired-2.fq.gz ... ", $geno, $r, $genotypes{$geno}->{"lib"});
			if (!exists $genotypes{$geno}->{"reps"}->{$r}->{"paired"}->{2}) {
				print STDERR sprintf("NOT PRESENT\n");
				$ok = false;
			} # end of if statement
			else {
				print STDERR sprintf("PRESENT\n");
			} # End of else statement

			print STDERR sprintf("        - %s.%s.%s.trimmed.repmasked-singletons-1.fq.gz ... ", $geno, $r, $genotypes{$geno}->{"lib"});
			if (!exists $genotypes{$geno}->{"reps"}->{$r}->{"singletons"}->{1}) {
				print STDERR sprintf("NOT PRESENT\n");
				$ok = false;
			} # end of if statement
			else {
				print STDERR sprintf("PRESENT\n");
			} # End of else statement
			
			print STDERR sprintf("        - %s.%s.%s.trimmed.repmasked-singletons-2.fq.gz ... ", $geno, $r, $genotypes{$geno}->{"lib"});
			if (!exists $genotypes{$geno}->{"reps"}->{$r}->{"singletons"}->{2}) {
				print STDERR sprintf("NOT PRESENT\n");
				$ok = false;
			} # end of if statement
			else {
				print STDERR sprintf("PRESENT\n");
			} # End of else statement
		} # end of foreach statement

		print STDERR sprintf("    + Generating assembly script '%s' ... ", $script_name);
		if (!$ok) {
			print STDERR sprintf("FAILED\n");
		} # end of if statement
		else {	# Generate script
			my $ofh = new FileHandle();
			open ($ofh, sprintf(">%s", $script_name)) or die("Cannot generate script\n");
			print $ofh sprintf("#!/bin/bash\n");
			print $ofh sprintf("\n");
			print $ofh sprintf("genotype=\"%s\"\n", $geno);
			print $ofh sprintf("output_dir=\"%s\"\n", $outputDir);
			print $ofh sprintf("proc=%s\n", $proc);
			print $ofh sprintf("aligner=\"%s\"\n", DEFAULT_ALIGNER);
			print $ofh sprintf("k=%s\n", $k);
			print $ofh sprintf("c=%s\n", DEFAULT_C);
			print $ofh sprintf("s=%s\n", DEFAULT_S);
			print $ofh sprintf("n=%s\n", DEFAULT_N);
			print $ofh sprintf("m=%s\n", DEFAULT_M);
			print $ofh sprintf("p=%s\n", DEFAULT_P);
			print $ofh sprintf("\n");
			print $ofh sprintf("# Decompressing files\n");

			my (@PE, @lblPE, @PEFiles, @SE, @lblSE, @SEFiles);
			my ($PECount, $SECount) = (0, 0);
			foreach my $r (@reps) {
				foreach my $t (sort {$a cmp $b} keys %{ $genotypes{$geno}->{"reps"}->{$r} }) {
					my $f1 = $genotypes{$geno}->{"reps"}->{$r}->{$t}->{1};
					$f1 =~ s/\.gz$//g;
					print $ofh sprintf("gzip -cd %s/%s.gz > %s\n", $inputDir, $f1, $f1);
					
					my $f2 = $genotypes{$geno}->{"reps"}->{$r}->{$t}->{2};
					$f2 =~ s/\.gz$//g;
					print $ofh sprintf("gzip -cd %s/%s.gz > %s\n", $inputDir, $f2, $f2);
						
					if ($t =~ m/^paired$/i) {	# Paired-end
						my $id = sprintf("PE%s", ++$PECount);
						push(@lblPE, $id);
						push(@PE, sprintf("%s=\"%s %s\"", $id, $f1, $f2));
						push(@PEFiles, $f1);
						push(@PEFiles, $f2);
					} # end of if statement
					elsif ($t =~ m/^singletons$/i) {	# Singletons
						my $id = sprintf("SE%s", ++$SECount);
						push(@lblSE, $id);
						push(@SE, sprintf("%s=\"%s %s\"", $id, $f1, $f2));
						push(@SEFiles, $f1);
						push(@SEFiles, $f2);
					} # end of else if statement
					else {
						die("Invalid Type: $t\n");
					} # end of else statement
				} # end of for each statement
			} # end of for each statement

			print $ofh sprintf("\n");
			print $ofh sprintf("# Paired-End libraries\n");
			print $ofh sprintf("%s\n", join("\n", @PE));
			
			print $ofh sprintf("\n");
			print $ofh sprintf("# Single-End libraries\n");
			print $ofh sprintf("%s\n", join("\n", @SE));

			print $ofh sprintf("\n");
			print $ofh sprintf("# Assembly\n");
			print $ofh sprintf("time -v -o \"\$genotype\".k\"\$k\"-\"\$aligner\".resources abyss-pe default clean \\\n");
			print $ofh sprintf("\taligner=\"\$aligner\" np=\"\$proc\" c=\"\$c\" s=\"\$s\" k=\"\$k\" n=\"\$n\" j=\"\$proc\" m=\"\$m\" p=\"\$p\" \\\n");
			print $ofh sprintf("\tname=\"\$genotype\".k\"\$k\"-\"\$aligner\" \\\n");
			print $ofh sprintf("\tlib=\"%s\" ", join(" ", @lblPE));
			for (my $i=0; $i < scalar(@lblPE); $i++) {
				print $ofh sprintf(" ") if ($i != 0);
				print $ofh sprintf("%s=\"\$%s\"", $lblPE[$i], $lblPE[$i]);
			} # end of for loop
			print $ofh sprintf(" \\\n");

			print $ofh sprintf("\tse=\"%s\" ", join(" ", map { sprintf("\$%s", $_); } @lblSE));
			print $ofh sprintf(" \\\n");
			
			print $ofh sprintf("\t>> \"\$genotype\".k\"\$k\"-\"\$aligner\".log 2> \"\$genotype\".k\"\$k\"-\"\$aligner\".err\n");

			print $ofh sprintf("\n");
			print $ofh sprintf("# Rename contigs\n");
			print $ofh sprintf("rename_abyss_assembly.pl -i \"\$genotype\".k\"\$k\"-\"\$aligner\"-contigs.fa -o \"\$genotype\".k\"\$k\"-\"\$aligner\"-contigs.gte300.fas\n");

			print $ofh sprintf("\n");
			print $ofh sprintf("# Organization\n");
			print $ofh sprintf("mkdir -p \"\$output_dir\"/\"\$genotype\"\n");
			print $ofh sprintf("mv %s \"\$output_dir\"/\"\$genotype\"\n", $script_name);
			print $ofh sprintf("mv \"\$genotype\".k\"\$k\"-\"\$aligner\".* \"\$output_dir\"/\"\$genotype\"\n");
			print $ofh sprintf("mv \"\$genotype\".k\"\$k\"-\"\$aligner\"-* \"\$output_dir\"/\"\$genotype\"\n");

			print $ofh sprintf("\n");
			print $ofh sprintf("# Cleanup\n");
			foreach my $pefile (@PEFiles) {
				print $ofh sprintf("rm -r %s\n", $pefile);
			} # end of foreach statement
			foreach my $sefile (@SEFiles) {
				print $ofh sprintf("rm -r %s\n", $sefile);
			} # end of foreach statement

			close ($ofh);
			print STDERR sprintf("DONE\n");
			
			print STDERR sprintf("    + Running assembly '%s' ... ", $script_name);
			my $command = sprintf("sh %s", $script_name);
			system($command);
			print STDERR sprintf("DONE\n");

			my $end_time = timelocal(localtime(time));
			print STDERR sprintf("    + End Time: %s\n", scalar(localtime(time)));
			print STDERR sprintf("    + Run-Time: %s\n", &formatTime($end_time - $start_time));
		} # End of else statement
	} # end of else statement
	print STDERR sprintf("\n");
} # end of for loop

print STDERR sprintf("\n");


