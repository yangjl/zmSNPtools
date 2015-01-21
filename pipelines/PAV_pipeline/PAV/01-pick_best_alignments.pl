#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Schnablelab::Tools;

use constant TAIL_CUTOFF => 5;

sub GetAttributes {
	my $str = $_[0];
	my ($name, $length, $ident, $cvg, $tail_f, $tail_e, $overhang, $aligned, $start, $end, $score);

	if ($str =~ m/Length: (\d+) bp/) {
		$length = $1;
	} # end of if statemetn

	if ($str =~ m/Identity: (\d+\.\d+)/) {
		$ident = $1;
	} # end of if statement

	if ($str =~ m/Coverage: (\d+\.\d+)/) {
		$cvg = $1;
	} # end of if statemetn

	if ($str =~ m/Tails Front: (\d+) bp/) {
		$tail_f = $1;
	} # end of if statement

	if ($str =~ m/Tails End: (\d+) bp/) {
		$tail_e = $1;
	} # end of if statement

	if ($str =~ m/Overhang: (\d+) bp/) {
		$overhang = $1;
	} # end of if statemetn

	if ($str =~ m/Aligned: (\d+) bp/) {
		$aligned = $1;
	} # end of if statement

	if ($str =~ m/Target=(\S+) (\d+) (\d+)/) {
		$name = $1;
		$start = $2;
		$end = $3;
	} # end of if statement

	if ($str =~ m/Score: (\S+);Target/) {
		$score = $1;
	} # end of if statemetn

	return ($name, $length, $ident, $cvg, $tail_f, $tail_e, $overhang, $aligned, $start, $end, $score);
} # end of sub GetAttributes

my (@files, $cutoff);
my $result = &GetOptions("files|f=s{1,}" => \@files, "tail|t:i{1}" => \$cutoff);

unless ($result && scalar(@files) > 0) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --files <*.gff3>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

$cutoff = TAIL_CUTOFF if (!defined($cutoff) || $cutoff !~ m/^\d+$/);

print STDERR sprintf("\n");
my %contigs;
foreach my $f (@files) {
	my $source;
	if ($f =~ m/ZmB73_RefGen_v3\.pieces\+mito\+chlo/) {
		$source = "AGPv3";
	} # end of if statemetn
	elsif ($f =~ m/ZmB73_RefGen_v2\.pieces\+mito\+chlo/) {
		$source = "AGPv2";
	} # end of else if statemetn
	elsif ($f =~ m/ZmB73_RefGen_v1\.pieces/) {
		$source = "AGPv1";
	} # end of else if statement
	elsif ($f =~ m/WUGSC_BACs-20130313\.pieces/) {
		$source = "WUGSC";
	} # end of else if statemetn
	elsif ($f =~ m/Nat_Genet_B73\.novel/) {
		$source = "Nat_Genet";
	} # End of else if statemetn
	elsif ($f =~ m/MAGIv3\.1\.contigs/) {
		$source = "MAGIv3.1";
	} # end of else if statement
	elsif ($f =~ m/MAGIv4\.0\.contigs/) {
		$source = "MAGIv4.0";
	} # end of else if statemetn
	else {
		$source = "OTHER";
	} # end of else statemetn

	if (defined($source)) {
		print STDERR sprintf("o Procesing '%s' ... ", $f);
		my $fh = new FileHandle();
		open ($fh, $f) or die("Cannot open file\n");
		while (<$fh>) {
			chomp;
			if (length($_) != 0 && $_ !~ m/^#/) {
				my @fields = split(/\t/, $_);
				my $ref = $fields[0];
				my $ref_start = &min($fields[3], $fields[4]);
				my $ref_end = &max($fields[3], $fields[4]);
				my ($name, $length, $ident, $cvg, $tail_f, $tail_e, $overhang, $aligned, $start, $end, $score) = &GetAttributes($fields[$#fields]);
				
				if (!exists $contigs{$name} || $score > $contigs{$name}->{"score"}) {
					$contigs{$name}->{"length"} = $length;
					$contigs{$name}->{"ref"} = $ref;
					$contigs{$name}->{"ref_start"} = $ref_start;
					$contigs{$name}->{"ref_end"} = $ref_end;
					$contigs{$name}->{"source"} = $source;
					$contigs{$name}->{"ident"} = $ident;
					$contigs{$name}->{"cvg"} = $cvg;
					$contigs{$name}->{"tail_f"} = $tail_f;
					$contigs{$name}->{"tail_e"} = $tail_e;
					$contigs{$name}->{"overhang"} = $overhang;
					$contigs{$name}->{"aligned"} = $aligned;
					$contigs{$name}->{"start"} = $start;
					$contigs{$name}->{"end"} = $end;
					$contigs{$name}->{"score"} = $score;
				} # end of if statement
			} # end of if statement
		} # end of while loop
		close ($fh);
		print STDERR sprintf("DONE\n");
	} # end of if statement
} # end of for each statement

my %categories;
my $total = 0;
printf("# NAME\tLENGTH\tALIGN_SRC\tREF\tREF_START\tREF_END\tIDENT\tCVG\tTAIL_F\tTAIL_E\tOVERHANG\tALIGNED\tALIGN_START\tALIGN_END\tBREAK_LEFT\tBREAK_RIGHT\tCATEGORY\tTOTAL_TAIL\tTAIL_PER\n");
foreach my $name (keys %contigs) {
	printf("%s", $name);
	printf("\t%s", $contigs{$name}->{"length"});
	printf("\t%s", $contigs{$name}->{"source"});
	printf("\t%s", $contigs{$name}->{"ref"});
	printf("\t%s", $contigs{$name}->{"ref_start"});
	printf("\t%s", $contigs{$name}->{"ref_end"});
	printf("\t%s", $contigs{$name}->{"ident"});
	printf("\t%s", $contigs{$name}->{"cvg"});
	printf("\t%s", $contigs{$name}->{"tail_f"});
	printf("\t%s", $contigs{$name}->{"tail_e"});
	printf("\t%s", $contigs{$name}->{"overhang"});
	printf("\t%s", $contigs{$name}->{"aligned"});
	printf("\t%s", $contigs{$name}->{"start"});
	printf("\t%s", $contigs{$name}->{"end"});
	
	if ($contigs{$name}->{"tail_f"} != 0 && $contigs{$name}->{"start"} - 1 != 0) {
		printf("\t%s", $contigs{$name}->{"start"} - 1);
	} # end of if statement
	else {
		printf("\t---");
	} # end of if statement

	if ($contigs{$name}->{"tail_e"} != 0 && $contigs{$name}->{"end"} != $contigs{$name}->{"length"}) {
		printf("\t%s", $contigs{$name}->{"end"} + 1);
	} # end of if statement
	else {
		printf("\t---");
	} # end of else statement

	my $cat;
	if ($contigs{$name}->{"tail_f"} == 0 && $contigs{$name}->{"tail_e"} == 0) {
		$cat = sprintf("NO_TAIL");
	} # end of if statemetn
	elsif ($contigs{$name}->{"tail_f"} != 0 && $contigs{$name}->{"tail_e"} == 0) {
		if ($contigs{$name}->{"tail_f"} <= $cutoff) {
			$cat = sprintf("FRONT/END_LESSTHAN_%s", $cutoff);
		} # end of if statemnet
		else {
			$cat = sprintf("FRONT/END_MORETHAN_%s", $cutoff);
		} # end of else statement
	} # end of else if statement
	elsif ($contigs{$name}->{"tail_f"} == 0 && $contigs{$name}->{"tail_e"} != 0) {
		if ($contigs{$name}->{"tail_e"} <= $cutoff) {
			$cat = sprintf("FRONT/END_LESSTHAN_%s", $cutoff);
		} # end of if statemnet
		else {
			$cat = sprintf("FRONT/END_MORETHAN_%s", $cutoff);
		} # end of else statement
	} # end of else if statemetn
	else {
		my $tail = $contigs{$name}->{"tail_f"} + $contigs{$name}->{"tail_e"};
		if ($tail <= $cutoff) {
			$cat = sprintf("BOTH_LESSTHAN_%s", $cutoff);
		} # endo f if statemetn
		else {
			$cat = sprintf("BOTH_MORETHAN_%s", $cutoff);
		} # end of else statemetn
	} # end of else statemetn
	
	printf("\t%s", $cat);
	if (exists $categories{$cat}) {
		$categories{$cat}++;
	} # end of if statemetn
	else {
		$categories{$cat} = 1;
	} # end of else statemetn
	$total++;

	my $tail = $contigs{$name}->{"tail_f"} + $contigs{$name}->{"tail_e"};
	printf("\t%s\t%s", $tail, $tail / $contigs{$name}->{"length"});

	printf("\n");
} # end of for each statement

my $notail = 0;
my $frontend = 0;
my $both = 0;

foreach my $c (sort {$a cmp $b} keys %categories) {
	printf("# %s\t%s\t%2.1f%%\n", $c, &formatNumber($categories{$c}), ($categories{$c} / $total) * 100);
	if ($c =~ m/^NO_TAIL$/i) {
		$notail += $categories{$c};
	} # End of if statemetn
	elsif ($c =~ m/^FRONT\/END/) {
		$frontend += $categories{$c};
	} # End of else if statemetn
	elsif ($c =~ m/^BOTH/) {
		$both += $categories{$c};
	} # End of else statement
} # end of foreach statemetn
printf("# TOTAL\t%s\n", &formatNumber($total));

printf("\n");
printf("# NO_TAIL: %s\n", &formatNumber($notail));
printf("# FRONT/END: %s\n", &formatNumber($frontend));
printf("# BOTH: %s\n", &formatNumber($both));

print STDERR sprintf("\n");
