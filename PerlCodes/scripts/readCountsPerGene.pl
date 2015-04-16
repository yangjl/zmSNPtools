#!/usr/bin/perl -w

# Author: Eddy Yeh
#
# Given parsed GSNAP tabular files, GSNAP gff3 files, and native files, get the read counts of each gene
#
# Change Log v1.9 (Eddy) - 2012.07.25
#  o Removed support for tabular alignment format
#  o Removed support for native alignment format (Internal native SNP calling files)
#
# Change Log v1.8 (Eddy) - 2012.07.20
#  o Added --fragment|--nofragment arguments. This argument controls the read count weight for paired-end
#    reads within a gene. If  the value is specified (--fragment) and when both mates of a fragment are within
#    the same gene, the tally count of the fragment is 1 instead of 2 (although both reads are within the same
#    gene). When --nofragment option is specified, it disables the fragment tally and allocates 2 read counts
#    for the entire fragment within the gene.
#
# Change Log v1.7 (Eddy) - 2012.04.23
#  o Fixed bug when attempting to count reads per gene using native files
#
# Change Log v1.6 (Eddy) - 2012.04.10
#  o Added support to count SNPs per gene from polymorphisms detected using our in-house SNP discovery
#    pileine (--snps)
# 
# Change Log v1.5 (Eddy) - 2011.12.16
#  o Fixed inconsistent chr0 chromosome present in Maize annotations to chrUNKNOWN before lookups
#  o Fixed bitwise operator in sam output flag field to capture read orientation in alignment more accurately
#
# Change Log v1.4 (Eddy) - 2011.08.27
#  o Fixed small bug causing some genes to have incorrect number of read counts due to
#    the way that the array of positions were sorted.
#  o Reads having coordinatates that spawn across a gene are no longer valid. For a read to be
#    considered inside a gene, either the start or end positions have to be present in the gene
#
# Change Log v1.3 (Eddy) - 2011.07.31
#  o Change implementation. Instead of splitting genes and reads into separate files, keep them in memory
#    for processing
#
# Change Log v1.2 (Eddy) - 2011.05.24
#  o Added support to perform read counts per gene having read alignments in SAM format
#
# Change Log v1.1 (Eddy) - 2011.05.12
#  o Added support of Sanzhen's gene table in addition
# 
# Change Log v1.0 (Eddy) - 2011.03.31
#  o Initial implementation of this script
#

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Time::Local;
use Storable qw(freeze thaw);
use POSIX qw(ceil floor);

use constant true => 1;
use constant false => 0;
use constant DEFAULT_MIN_READ_COUNT => 0;
use constant DEFAULT_FRAGMENT_COUNT => 1;	# Count to be tallied towards paired-end reads within the same gene
use constant LIMIT => 1000;
use constant VERSION => "1.9 (2012.07.25)";

$| = true;

sub min {
    my @sorted = sort {$a <=> $b} @_; 
    return shift(@sorted);
} # End of sub min

sub max {
    my @sorted = sort {$a <=> $b} @_; 
    return pop(@sorted);
} # End of sub max

sub trim {
    my $str = $_[0];
    
    if (defined($str)) {
        $str =~ s/^\s+//g;
        $str =~ s/^\t+//g;
        $str =~ s/\s+$//g;
        $str =~ s/\t+$//g;
    } # End of if statement
    else {
        $str = "";
    } # End of else statemnet

    return $str;
} # End of sub trim

sub removeComma {
	my $str = $_[0];
	
	$str =~ s/,//g;
	return $str;
} # end of sub removeComma

sub addPadding {
	my ($str, $length) = @_;
	my $padded = $str;

	for (my $i=length($str); $i < $length; $i++) {
		$padded .= " ";
	} # end of for loop

	return $padded;
} # end of sub addPadding

sub formatNumber {
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;

    return $_;
} # End of sub formatNumber

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

sub fragmentOverlap {
    my ($istart, $iend, $jstart, $jend) = @_;
    my $returnVal = false;  # False default

    if (($jstart >= $istart && $jstart <= $iend) ||
        ($istart >= $jstart && $istart <= $jend)) {
        $returnVal = true;
    } # End of if statement

    return $returnVal;
} # End of sub fragmentOverlap

sub getName {
	my @tokens = split(/;/, $_[0]);
	my $name;
	my $found = false;
	for (my $i=0; $i < scalar(@tokens) && !$found; $i++) {
		if ($tokens[$i] =~ m/^Name=/i) {
			$name = $';
			$found = true;
		} # End of if statement
	} # End of for loop

	return $name;
} # End of getName

# Find the index in $positions in which the read overlaps with the given gene coordinates
sub findIndex {
	my ($gene_start, $gene_end, $positions) = @_;
	my $start = 0;		# start index
	my $end = scalar(@{ $positions }) - 1;	# end index
	my $index;
	my ($r_coord1, $r_coord2, $r_name);

	while ($start <= $end && !defined($index)) {
		# Check if value of start index is enclosed inside gene
		($r_coord1, $r_coord2, $r_name) = split(/;/, $positions->[$start]);

		if ($r_coord1 >= $gene_start && $r_coord1 <= $gene_end) {	# found in start index
			$index = $start;
		} # end of if statement
		else {
			($r_coord1, $r_coord2, $r_name) = split(/;/, $positions->[$end]);
		
			# Check if value of end index is enclosed inside gene
			if ($r_coord1 >= $gene_start && $r_coord1 <= $gene_end) {	# found in end index
				$index = $end;
			} # end of if statement
			else {
				# Check if value of middle is enclosed inside gene
				my $middle = int(($end + $start) / 2);
				($r_coord1, $r_coord2, $r_name) = split(/;/, $positions->[$middle]);

				if ($r_coord1 >= $gene_start && $r_coord1 <= $gene_end) {	# found in middle index
					$index = $middle;
				} # end of if statement
				else {
					if ($r_coord1 < $gene_start) {
						# middle pointer points to reads that are to the left of the gene coordinates
						# so, we need to discard the left half and keep searching on the right side
						$start = $middle + 1;
					} # end of if statement
					elsif ($r_coord1 > $gene_end) {
						# middle pointer points to reads that are to the right of the gene coordinates
						# so, we need to discard right half and keep searching on the left side
						$end = $middle - 1;
					} # end of elsif statement
					else {
						die("Cannot discard either left or right half of positions\n");
					} # End of else statement
				} # End of if statement
			} # End of else statement
		} # End of else statemnet
	} # End of while loop

	if (defined($index)) {	# Search successful, move all the way to the left
		my $stop = false;
		for (my $i=$index - 1; $i >= 0 && !$stop; $i--) {
			($r_coord1, $r_coord2, $r_name) = split(/;/, $positions->[$i]);
			if ($r_coord1 >= $gene_start && $r_coord1 <= $gene_end) {
				$index--;	# Move to the left
			} # end of if statmenet
			else {
				$stop = true;
			} # end of else statement
		} # end of for loop
	} # end of if statmenet

	return $index;
} # End of sub findIndex

# Recompute read counts treating each paired-end read as 1 instead of 2. Basically the subroutine below
# removes duplicated read IDs from the array of names passed to the subroutine
sub recomputeReadCounts {
	my @names = @_;
    my %seen;

    foreach my $n (@names) {
        if (exists $seen{$n}) {
            $seen{$n}++;
        } # end of if statement
        else {
            $seen{$n} = 1;
        } # end of else statement
    } # end of for each statement

    return scalar(keys %seen);
} # end of sub recomputeReadCounts

# Count the number of genes from the genes array reference
sub countReads {
	my ($ofh, $genes, $array_start, $array_end, $num, $fragCount, $showNames, $printStats, $read_start_positions, $read_end_positions) = @_;
	my @read_names;
	my %ids;	# To keep track of unique read IDs
	my %master;	# To keep track of all the read IDs in gene space

	for (my $i=$array_start; $i <= $array_end; $i++) {
		@read_names = ();	# Clear read names
		%ids = ();			# Clear Ids
		my ($gid, $c, $gene_start, $gene_end, $strand) = split(/\t/, $genes->[$i]);

		# Find the index from read_start_positions that overlap with the gene coordinates
		my $index = &findIndex($gene_start, $gene_end, $read_start_positions);

		# If index is defined, that means $read_start_positions->[$index] is inside
		# the gene; Otherwise search failed (no reads inside gene)
		if (defined($index)) {
			my $ok = true;
			my ($r_coord1, $r_coord2, $r_name);

			for (my $i=$index; $i < scalar(@{ $read_start_positions }) && $ok; $i++) {	# Keep adding acceptable reads
				($r_coord1, $r_coord2, $r_name) = split(/;/, $read_start_positions->[$i]);
				if ($r_coord1 >= $gene_start && $r_coord1 <= $gene_end) {
					my $r_start = &min($r_coord1, $r_coord2);
					my $r_end = &max($r_coord1, $r_coord2);
					$ids{$r_start}->{$r_end}->{$r_name} = true;
				} #e nd of if statement
				else {
					$ok = false;
				} # end of else statmeent
			} # end of for loop
		} # End of if statement

		# Find the index from read_end_positions that overlap with the gene coordinates
		$index = &findIndex($gene_start, $gene_end, $read_end_positions);

		# If index is defined, that means $read_end_positions->[$index] is inside
		# the gene; Otherwise search failed (no reads inside gene)
		if (defined($index)) {
			my $ok = true;
			my ($r_coord1, $r_coord2, $r_name);

			for (my $i=$index; $i < scalar(@{ $read_end_positions }) && $ok; $i++) {	# Keep adding acceptable reads
				($r_coord1, $r_coord2, $r_name) = split(/;/, $read_end_positions->[$i]);
				if ($r_coord1 >= $gene_start && $r_coord1 <= $gene_end) {
					my $r_start = &min($r_coord1, $r_coord2);
					my $r_end = &max($r_coord1, $r_coord2);
					$ids{$r_start}->{$r_end}->{$r_name} = true;
				} # end of if statement
				else {
					$ok = false;
				} # end of else statement
			} # end of for loop
		} # end of if statement

		# Calculating read names
		foreach my $s (sort {$a <=> $b} keys %ids) {
			foreach my $e (sort {$a <=> $b} keys %{ $ids{$s} }) {
				foreach my $n (sort {$a cmp $b} keys %{ $ids{$s}->{$e} }) {
					push(@read_names, $n);

					# To keep track of read IDs in gene space
					$master{sprintf("%s_%s-%s", $n, $s, $e)} = true;
				} # end of for each statement
			} # End of for each statement
		} # End of for each statemnet

		if (scalar(@read_names) == 0 && $num == 0) {
			if ($showNames) {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				                   $gid, $c, $gene_start, $gene_end, $strand, 0, "-");
			} # end of if statement
			else {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",
				                   $gid, $c, $gene_start, $gene_end, $strand, 0);
			} # End of else statement
		} # end of if statement
		elsif (scalar(@read_names) >= $num) {
			# Read names are sorted based on their start/end coordinates
			my $readCounts;

            if ($fragCount) {
                # Need to recompute read counts treating each paired-end fragment as 1 instead of 2
                $readCounts = &recomputeReadCounts(@read_names);
            } # End of if statement
            else {
                $readCounts = scalar(@read_names);
            } # end of else statement

			if ($showNames) {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				                   $gid, $c, $gene_start, $gene_end, $strand, $readCounts, join(";", @read_names));
			} # end of if statement
			else {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",
				                   $gid, $c, $gene_start, $gene_end, $strand, $readCounts);
			} # End of else statement
		} # End of if statemetn
		elsif (!defined($index) && $num == 0) {
			if ($showNames) {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				                   $gid, $c, $gene_start, $gene_end, $strand, 0, "-");
			} # end of if statement
			else {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",
				                   $gid, $c, $gene_start, $gene_end, $strand, 0);
			} # End of else statement
		} # end of if statement
	} # end of for each statement

	my $reads_in_genes = keys %master;
	
	if ($printStats) {
		print $ofh sprintf("# TOTAL READS IN GENES : %s\n", &formatNumber($reads_in_genes));
	} # end of if statement

	return $reads_in_genes;
} # end of sub countReads

# Format the progress string
sub printProgress {
    my ($old_progress, $new_progress) = @_; 
    my $end_of_line = false;

    if ($new_progress =~ m/\n$/) {
        $end_of_line = true;
        chomp($new_progress);
    } # End of if statement

    if (defined($old_progress)) {
        for (my $i=0; $i < length($old_progress); $i++) {
            print STDERR sprintf("\b"); # Clear previous progress
        } # End of for loop
    } # End of if statemnet

    print STDERR sprintf("%s", $new_progress);

    # adjusting/padding text with white-space
    if (defined($old_progress) && length($new_progress) < length($old_progress)) {
        for (my $i=length($new_progress); $i < length($old_progress); $i++) {
            print STDERR sprintf(" ");
        } # End of for loop

        for (my $i=length($new_progress); $i < length($old_progress); $i++) {
            print STDERR sprintf("\b");
        } # End of for loop
    } # End of if statement

    if ($end_of_line) {
        print STDERR sprintf("\n");
        $new_progress .= "\n";
    } # End of if statement

    return $new_progress;
} # end of sub printProgress

# --------------------- MAIN PROGRAM STARTS HERE ---------------------
my ($genes_file, $sl_genes_file, @sam_files, @gff3_files, @snp_files, $output_file, $sample_name, $num, $showNames);
my $fragCount;

my $result = &GetOptions("genes|g:s{1}" => \$genes_file,
                         "gsl:s{1}" => \$sl_genes_file,
						 "gff3:s{1,}" => \@gff3_files,
						 "sam:s{1,}" => \@sam_files,
						 "snps:s{1}" => \@snp_files,
						 "output|out|o=s{1}" => \$output_file,
						 "sample|s:s{1}" => \$sample_name,
						 "num|n:i{1}" => \$num,
						 "fragment|f!" => \$fragCount,
						 "names!" => \$showNames);

unless ($result && (defined($genes_file) || defined($sl_genes_file))  && defined($output_file) && defined($sample_name)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("*** INCORRECT NUMBER OF ARGUMENTS ***\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("   perl %s [--genes|-g <genes.gff3>|--gsl <sanzhen's genes table>] --sample|-s <sample name> \\\n", $0);
	print STDERR sprintf("               --output|-out|-o <output file> [OPTIONS]\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --genes|-g <genes.gff3>                : Path to genes gff3 file with 'gene' feature\n");
	print STDERR sprintf("   --gsl <sanzhen's genes table>          : Path to Sanzhen's genes table file\n");
	print STDERR sprintf("   --sample|-s <sample name>              : Sample name to be included in the header\n");
	print STDERR sprintf("   --output|-out|-o <output file>         : Path to the output file were counts will be saved\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS\n");
	print STDERR sprintf("   --gff3 <gff3 files>                    : Parsed read alignment files in gff3 format\n");
	print STDERR sprintf("   --sam <sam files>                      : Parsed read alignment files in SAM format\n");
	print STDERR sprintf("   --snps <snp files>                     : Putative SNPs/INDELs tabular output generated by in-house\n");
	print STDERR sprintf("                                            SNP discovery pipeline\n");
	print STDERR sprintf("   --names|--nonames                      : Enable/Disable printing of read names [DEFAULT: --nonames]\n");
	print STDERR sprintf("   --num <int>                            : Only display genes having at least this number\n");
	print STDERR sprintf("                                            of reads and/or SNPs/INDELs [DEFAULT: %s]\n", DEFAULT_MIN_READ_COUNT);
	print STDERR sprintf("   --fragment|--nofragment                : This argument controls the read count weight for paired-end\n");
	print STDERR sprintf("                                            eads within a gene. If  the value is specified (--fragment) and \n");
	print STDERR sprintf("                                            when both mates of a fragment are within the same gene, the tally\n");
	print STDERR sprintf("                                            count of the fragment is 1 instead of 2 (although both reads are\n");
	print STDERR sprintf("                                            within the same gene). When --nofragment option is specified, it\n");
	print STDERR sprintf("                                            disables the fragment tally and allocates 2 read counts for the\n");
	print STDERR sprintf("                                            entire fragment within the gene [DEFAULT: --fragment]\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("VERSION: %s\n", VERSION);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

# Check if alignment files were provided
if (scalar(@sam_files) + scalar(@snp_files) + scalar(@gff3_files) == 0) {
	print STDERR sprintf("\n");
	print STDERR sprintf("ERROR: No files were specified containing alignment reads. Please see\n");
	print STDERR sprintf("       --sam, --gff3 and/or --snps options\n");
	print STDERR sprintf("\n");
	exit();
} # End of if statement

# Setting default parameters if missing
$num = DEFAULT_MIN_READ_COUNT if (!defined($num) || $num !~ m/^\d+$/);
$fragCount = true if (!defined($fragCount));
$showNames = false if (!defined($showNames));

my $start_time = timelocal(localtime(time));
my $counter = 0;
my $progress = "Please Wait ...";
my ($fh, $ofh);
my $total_reads = 0;
my $total_genes = 0;

print STDERR sprintf("\n");

my %master_genes;
my %master_reads;
my %master_counts;

# Processing standard gff3 file containing gene coordinates
if (defined($genes_file)) {
	$counter = 0;
	$progress = "Please Wait ...";
	print STDERR sprintf("  o Sorting genes file (%s) by chromosome : %s", $genes_file, $progress);
	$fh = new FileHandle();
	open ($fh, $genes_file) or die("Cannot open genes file\n");

	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/#/) {
			my ($ref, $source, $feature, $pos1, $pos2, $score, $strand, $frame, $group) = split(/\t/, $_);
			
			if ($feature =~ m/^gene$/i) {
				my ($chr, $start, $end, $name);
				
				# Fixing chr
				if ($ref =~ m/^unknown$/i) {
					$chr = "chrUNKNOWN";
				} # end of if statemetn
				elsif ($ref =~ m/^(\d+)$/ || $ref =~ m/^(Mt)$/i || $ref =~ m/^(Pt)$/i) {
					$chr = sprintf("chr%s", $1);
				} # end of if statemnet
				else {
					$chr = $ref;
				} # end of else statement

				$start = &min($pos1, $pos2);
				$end = &max($pos1, $pos2);
				$name = &getName($group);

				if (exists $master_genes{$chr}) {
					$master_genes{$chr} .= sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
					$master_counts{$chr}->{"gene_counts"}++;
				} # end of if statement
				else {
					$master_genes{$chr} = sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
					$master_counts{$chr}->{"gene_counts"} = 1;
				} # end of if statement

				$total_genes++;
			
				if (++$counter % LIMIT == 0) {
					$progress = &printProgress($progress, sprintf("%s genes processed", &formatNumber($counter)));
				} # End of if statement
			} # end of if statement
		} # end of if statement
	} # End of while loop
	close ($fh);

	&printProgress($progress, sprintf("%s genes total\n", &formatNumber($counter)));
} # end of if statement

# Processing Sanzhen's genes table file
if (defined($sl_genes_file)) {
	$counter = 0;
	$progress = "Please Wait ...";
	print STDERR sprintf("  o Sorting genes file (%s) by chromosome : %s", $sl_genes_file, $progress);
	
	$fh = new FileHandle();
	open ($fh, $sl_genes_file) or die("Cannot open Sanzhen's genes file\n");

	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^GeneID\tRef\tChr/) {
			my ($name, $ref, $chr, $strand, $pos1, $pos2, $exon) = split(/\t/, $_);
			my $start = &min($pos1, $pos2);
			my $end = &max($pos1, $pos2);

			# Fixing chr
			if ($chr =~ m/^unknown$/i || $chr =~ m/^chr0$/i || $chr =~ m/^0$/i) {
				$chr = "chrUNKNOWN";
			} # end of if statemetn
			elsif ($chr =~ m/^(\d+)$/ || $chr =~ m/^(Mt)$/i || $chr =~ m/^(Pt)$/i) {
				$chr = sprintf("chr%s", $1);
			} # end of if statemnet

			if (exists $master_genes{$chr}) {
				$master_genes{$chr} .= sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
				$master_counts{$chr}->{"gene_counts"}++;
			} # end of if statement
			else {
				$master_genes{$chr} = sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
				$master_counts{$chr}->{"gene_counts"} = 1;
			} # end of if statement

			$total_genes++;
		
			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s genes processed", &formatNumber($counter)));
			} # End of if statement
		} # end of if statement
	} # end of while loop
	close ($fh);
	&printProgress($progress, sprintf("%s genes total\n", &formatNumber($counter)));
} # end of if statement

print STDERR sprintf("  o Sorting read and polymorphism features by chromosome\n");

# Processing reads in gff3 format
foreach my $f (@gff3_files) {
	$counter = 0;
	$progress = "Please Wait ...";
	print STDERR sprintf("     + %s : %s", $f, $progress);
	$fh = new FileHandle();
	open ($fh, $f) or die("Cannot open gff3 file '$f' for reading\n");

	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my ($ref, $source, $feature, $pos1, $pos2, $score, $strand, $frame, $group) = split(/\t/, $_);
			
			if ($feature !~ m/^(fragment|HSP)$/i) {
				my ($chr, $start, $end, $name);

				# Fixing chr
				if ($ref =~ m/^unknown$/i || $ref =~ m/^chr0$/i || $ref =~ m/^0$/i) {
					$chr = "chrUNKNOWN";
				} # end of if statemetn
				elsif ($ref =~ m/^(\d+)$/) {
					$chr = sprintf("chr%d", $1);
				} # end of if statemnet
				else {
					$chr = $ref;
				} # End of else statement

				$start = &min($pos1, $pos2);
				$end = &max($pos1, $pos2);
				$name = &getName($group);

				if (exists $master_reads{$chr}) {
					$master_reads{$chr} .= sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
					$master_counts{$chr}->{"read_counts"}++;
				} # end of if statemnet
				else {
					$master_reads{$chr} = sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
					$master_counts{$chr}->{"read_counts"} = 1;
				} # end of else statement

				$total_reads++;

				if (++$counter % LIMIT == 0) {
					$progress = &printProgress($progress, sprintf("%s entries processed", &formatNumber($counter)));
				} # end of if statement
			} # End of if statement
		} # End of if statement
	} # end of while loop
	close ($fh);

	&printProgress($progress, sprintf("%s reads total\n", &formatNumber($counter)));
} # End of for each statement

# Processing reads in SAM format
foreach my $f (@sam_files) {
	$counter = 0;
	$progress = "Please Wait ...";
	print STDERR sprintf("     + %s : %s", $f, $progress);
	$fh = new FileHandle();
	open ($fh, $f) or die("Cannot open SAM file '$f' for reading\n");

	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^@/) {
			my ($name, $flag, $ref, $pos, $score, $cigar) = split(/\t/, $_);
			my ($chr, $start, $end, $strand);

			# Fixing chr
			if ($ref =~ m/^unknown$/i || $ref =~ m/^chr0/i || $ref =~ m/^0$/i) {
				$chr = "chrUNKNOWN";
			} # end of if statemetn
			elsif ($ref =~ m/^(\d+)$/) {
				$chr = sprintf("chr%d", $1);
			} # end of if statemnet
			else {
				$chr = $ref;
			} # End of else statement

			# Bitwise and 0x10 flag to determine strand
			# Operation below returns >0 if 0x10 present and 0 if 0x10 not present
			if ($flag & 0x10) {
				$strand = "-";
			} # end of if statement
			else {
				$strand = "+";
			} # end of else statement

			$start = $pos;

			# Compute end coordinate from cigar string
			$end = $start;
			foreach my $s (split(/(\d+[MIDNSHP=X])/, $cigar)) {
				if ($s =~ m/^(\d+)([M|D|N|=|X])$/i) {
					$end += $1;
				} # End of if statement
			} # end of for each statement
			$end -= 1 if ($start != $end);

			if (exists $master_reads{$chr}) {
				$master_reads{$chr} .= sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
				$master_counts{$chr}->{"read_counts"}++;
			} # end of if statemnet
			else {
				$master_reads{$chr} = sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
				$master_counts{$chr}->{"read_counts"} = 1;
			} # end of else statement

			$total_reads++;

			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s reads processed", &formatNumber($counter)));
			} # End of if statement
		} # End of if statement
	} # End of while loop
	close ($fh);

	&printProgress($progress, sprintf("%s reads total\n", &formatNumber($counter)));
} # end of for each statement

# Processing SNPs in tabular form generated by our in-house SNP calling pipeline
foreach my $f (@snp_files) {
	$counter = 0;
	$progress = "Please Wait ...";
	print STDERR sprintf("     + %s : %s", $f, $progress);
	$fh = new FileHandle();
	open ($fh, $f) or die("Cannot open SNP/INDEL file '$f' for reading\n");

	while (<$fh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my ($ref, $pos, $zygozity, $type, $ref_allele, $alt_allele, @junk) = split(/\t/, $_);
			my ($chr, $name, $start, $end, $strand);

			# Fixing chr
			if ($ref =~ m/^unknown$/i || $ref =~ m/^chr0/i || $ref =~ m/^0$/i) {
				$chr = "chrUNKNOWN";
			} # end of if statemetn
			elsif ($ref =~ m/^(\d+)$/) {
				$chr = sprintf("chr%d", $1);
			} # end of if statemnet
			else {
				$chr = $ref;
			} # End of else statement

			# Setting name, strand, start, and end positions
			$name = sprintf("SNP_%s_%s:%s/%s", $chr, $pos, $ref_allele, $alt_allele);
			$start = $pos;
			$end = $pos;
			$strand = "+";	# Useless, assuming strand is always 5' -> 3'

			if (exists $master_reads{$chr}) {
				$master_reads{$chr} .= sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
				$master_counts{$chr}->{"read_counts"}++;
			} # end of if statemnet
			else {
				$master_reads{$chr} = sprintf("%s\t%s\t%s\t%s\t%s\n", $name, $chr, $start, $end, $strand);
				$master_counts{$chr}->{"read_counts"} = 1;
			} # end of else statement

			$total_reads++;

			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s SNPs/INDELs processed", &formatNumber($counter)));
			} # End of if statement
		} # End of if statement
	} # End of while loop
	close ($fh);

	&printProgress($progress, sprintf("%s SNPs/INDELs total\n", &formatNumber($counter)));
} # end of for each statement

# Creating output file
$ofh = new FileHandle();
open ($ofh, sprintf(">%s", $output_file)) or die("Cannot open output file for writing\n");

# Print header information
print $ofh sprintf("# %s\n", scalar(localtime(time)));
print $ofh sprintf("#\n");
print $ofh sprintf("# Script: %s\n", $0);
print $ofh sprintf("# Version: %s\n", VERSION);
print $ofh sprintf("#\n");
print $ofh sprintf("# Genes File: %s\n", $genes_file) if (defined($genes_file));
print $ofh sprintf("# Genes File: %s\n", $sl_genes_file) if (defined($sl_genes_file));
print $ofh sprintf("# Sample Name: %s\n", $sample_name);
print $ofh sprintf("# Minimum number of reads per gene to show: %s\n", $num);
print $ofh sprintf("# Treat paired-end reads (fragment) as single-count: %s\n", $fragCount ? 
	"Yes (Paired-end reads within the same gene are counted as 1)" : "No (Paired-end reads within the same gene are counted as 2)");
print $ofh sprintf("# Show Read Names: %s\n", $showNames ? "Yes" : "No");

if (scalar(@gff3_files) > 0) {
	print $ofh sprintf("#\n");
	print $ofh sprintf("# GFF3 Files (%s):\n", &formatNumber(scalar(@gff3_files)));
	foreach my $f (@gff3_files) {
		print $ofh sprintf("#   o %s\n", $f);
	} # end of for each statement
} # End of if statement

if (scalar(@sam_files) > 0) {
	print $ofh sprintf("#\n");
	print $ofh sprintf("# SAM Files (%s):\n", &formatNumber(scalar(@sam_files)));
	foreach my $f (@sam_files) {
		print $ofh sprintf("#   o %s\n", $f);
	} # end of for each statement
} # End of if statement

if (scalar(@snp_files) > 0) {
	print $ofh sprintf("#\n");
	print $ofh sprintf("# SNP/INDEL Files (%s):\n", &formatNumber(scalar(@snp_files)));
	foreach my $f (@snp_files) {
		print $ofh sprintf("#   o %s\n", $f);
	} # end of for each statement
} # End of if statement

print $ofh sprintf("#\n");

# print header information
print $ofh sprintf("# Column Headers:\n");
print $ofh sprintf("# Gene ID\tChr.\tGene Start\tGene End\tStrand\tRead Counts");

if ($showNames) {
	print $ofh sprintf("\tRead IDs\n");
} # end of if statement
else {
	print $ofh sprintf("\n");
} # end of else statemnet

print $ofh sprintf("\n");

my $read_start_positions;	# Place holder for uniquely aligned read start positions
my $read_end_positions;		# Place holder for uniquely aligned read end positions
my $reads_in_genes = 0;		# Total reads in gene space
my $genes;					# Variable to contain gene positions in the chromosome
my @start_positions;		# Variable to contain sorted start positions of reads per chromosome  in ascending order
my @end_positions;			# Variable to contain sorted end positions of read per chromosome in ascending order
my @chromosomes = sort {$a cmp $b} keys %master_genes;

# Reads counting in genes
print STDERR sprintf("  o Reads counting in genes (%s reference sequences)\n", scalar(@chromosomes));

my $digits = length(sprintf("%s", &formatNumber(scalar(@chromosomes))));
my $c_num = 0;
foreach my $c (@chromosomes) {
	my $s_time = timelocal(localtime(time)); 
	undef($read_start_positions);		# Important, undefine variable to free up resources
	undef($read_end_positions);			# Important, undefine variable to free up resources
	undef($genes);						# Important, undefine variable to free up resources
	@start_positions = ();				# Important, reset read start positions of this chromosome
	@end_positions = ();				# Important, reset read end positions of this chomosome
	$counter = 0;
	$progress = "Please Wait ...";

	my $chr_genes = exists $master_counts{$c}->{"gene_counts"} ? $master_counts{$c}->{"gene_counts"} : 0;
	my $chr_reads = exists $master_counts{$c}->{"read_counts"} ? $master_counts{$c}->{"read_counts"} : 0;

	print STDERR sprintf("     + [ %\Q$digits\Es / %\Q$digits\Es ] Counting in '%s' [ %s genes ; %s features ]", 
	                     &formatNumber(++$c_num), &formatNumber(scalar(@chromosomes)), 
	                     $c, &formatNumber($chr_genes), &formatNumber($chr_reads));

	if ($chr_genes != 0 && $chr_reads != 0) {
		my $s_time = timelocal(localtime(time));
		print STDERR sprintf("\n");
		print STDERR sprintf("          - Loading read and polymorphism features : %s", $progress);
	
		my @lines = split(/\n/, $master_reads{$c});
		for (my $ii=0; $ii < scalar(@lines); $ii++) {
			my ($n, $c, $s, $e, $o) = split(/\t/, $lines[$ii]);
			
			# accumulate read start positions
			if (exists $read_start_positions->{$s}) {
				$read_start_positions->{$s} .= sprintf("\t%s;%s", $e, $n);
			} # end of if statement
			else {
				$read_start_positions->{$s} = sprintf("%s;%s", $e, $n);
			} # end of else statment

			# accumulate read end positions
			if (exists $read_end_positions->{$e}) {
				$read_end_positions->{$e} .= sprintf("\t%s;%s", $s, $n);
			} # end of if statement
			else {
				$read_end_positions->{$e} = sprintf("%s;%s", $s, $n);
			} # end of else statmenet

			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s (%2.1f%%) features processed", &formatNumber($counter), ($counter / $chr_reads) * 100));
			} # end of if statement
		} # end of for loop

		&printProgress($progress, sprintf("%s features loaded\n", &formatNumber($counter)));

		# Read genes file
		$counter = 0;
		$progress = "Please Wait ...";
		print STDERR sprintf("          - Loading genes : %s", $progress);

		undef($genes);
		@lines = split(/\n/, $master_genes{$c});
		for (my $ii=0; $ii < scalar(@lines); $ii++) {
			$genes->[$ii] = $lines[$ii];

			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s (%2.1f%%) genes processed", &formatNumber($counter), ($counter / $chr_genes) * 100));
			} # End of if statement
		} # end of for loop

		&printProgress($progress, sprintf("%s genes loaded\n", &formatNumber($counter)));

		# Preparing read start positions
		$counter = 0;
		$progress = "Please Wait ...";
		print STDERR sprintf("          - Preparing feature 'start' positions for lookup : %s", $progress);
		foreach my $s (sort {$a <=> $b} keys %{ $read_start_positions }) {
			foreach my $combo (split(/\t/, $read_start_positions->{$s})) {
				push(@start_positions, join(";", $s, $combo));
			} # end of for each statement

			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s (%2.1f%%) processed", &formatNumber($counter), ($counter / $chr_reads) * 100));
			} # end of if statement
		} # end of foreach statemnt
		&printProgress($progress, "COMPLETE\n");
		
		# Preparing read end positions
		$counter = 0;
		$progress = "Please Wait ...";
		print STDERR sprintf("          - Preparing feature 'end' positions for lookup : %s", $progress);
		foreach my $e (sort {$a <=> $b} keys %{ $read_end_positions }) {
			foreach my $combo (split(/\t/, $read_end_positions->{$e})) {
				push(@end_positions, join(";", $e, $combo));
			} # end of for each statement

			if (++$counter % LIMIT == 0) {
				$progress = &printProgress($progress, sprintf("%s (%2.1f%%) processed", &formatNumber($counter), ($counter / $chr_reads) * 100));
			} # end of if statement
		} # end of foreach statemnt
		&printProgress($progress, "COMPLETE\n");

		my $progress = "Please Wait ...";
		print STDERR sprintf("          - Processing Lookups : %s", $progress);

		$reads_in_genes += &countReads($ofh, $genes, 0, scalar(@{ $genes }) - 1, $num, $fragCount, $showNames, false, \@start_positions, \@end_positions);

		my $e_time = timelocal(localtime(time));
		&printProgress($progress, sprintf("%s\n", &formatTime($e_time - $s_time)));
	} # End of if statemnet
	else {
		if ($master_counts{$c}->{"gene_counts"} != 0 && (!exists $master_counts{$c}->{"read_counts"} || $master_counts{$c}->{"read_counts"} == 0)) {
			$counter = 0;
			$progress = "Please Wait ...";

			my @lines = split(/\n/, $master_genes{$c});

			print STDERR sprintf(" : %s", $progress);

			if ($num == 0) {
				for (my $ii=0; $ii < scalar(@lines); $ii++) {
					if ($showNames) {
						print $ofh sprintf("%s\t0\t-\n", $lines[$ii]);
					} # end of if statemnet
					else {
						print $ofh sprintf("%s\t0\n", $lines[$ii]);
					} # end of else statement

					if (++$counter % LIMIT == 0) {
						$progress = &printProgress($progress, sprintf("%s (%2.1f%%) genes processed", &formatNumber($counter), ($counter / $chr_genes) * 100));
					} # End of if statement
				} # End of for loop
			} # end of if statemnet

			&printProgress($progress, "NO FEATURES AVAILABLE\n");
		} # End of if statement
		else {
			print STDERR sprintf(" : IGNORED\n");
		} # end of else statement
	} # End of else statement
} # End of for each statement

if (defined($ofh)) {
	print $ofh sprintf("\n");
	print $ofh sprintf("# %s : %s\n", &addPadding("TOTAL GENES", 36), &formatNumber($total_genes));
	print $ofh sprintf("# %s : %s\n", &addPadding("TOTAL READS AND SNPs", 36), &formatNumber($total_reads));
	print $ofh sprintf("# %s : %s (%2.1f%%)\n", &addPadding("TOTAL READS AND SNPs/INDELs IN GENES", 36), &formatNumber($reads_in_genes), ($reads_in_genes / $total_reads) * 100); 

	print STDERR sprintf("  o %s : %s\n", &addPadding("TOTAL GENES", 36), &formatNumber($total_genes));
	print STDERR sprintf("  o %s : %s\n", &addPadding("TOTAL READS AND SNPs/INDELs", 36), &formatNumber($total_reads));
	print STDERR sprintf("  o %s : %s (%2.1f%%)\n", &addPadding("TOTAL READS AND SNPs/INDELs IN GENES", 36), &formatNumber($reads_in_genes), ($reads_in_genes / $total_reads) * 100);
} # end of if statement

my $end_time = timelocal(localtime(time));
print $ofh sprintf("# %s : %s\n", &addPadding("RUN TIME", 36), &formatTime($end_time - $start_time));
print STDERR sprintf("  o %s : %s\n", &addPadding("RUN TIME", 36), &formatTime($end_time - $start_time));
print STDERR sprintf("\n");
close ($ofh) if (defined($ofh));		# Close final output handler 
