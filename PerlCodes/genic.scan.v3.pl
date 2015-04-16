#!/usr/bin/perl -w

# genic.scan.pl

use strict;
use warnings;

use constant BINSIZE => 100000;

my (@line, $chr, %chr, $count, $bin, %count, $genecount);
# obtain genotyping info:
my $geno;
if ($ARGV[0] =~ /([\d\w]+)to/) {
	$geno = $1;
} else {
	$geno = $ARGV[0];
}

# load all reads info into hash
open(IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	@line = split(/\t/, $_);
	if (/chr(\d+)/) {
		$chr = $1;
		$count{$chr.$line[3]}++;
		$count = $count{$chr.$line[3]};
		$bin = int($line[3]/BINSIZE);
		$chr{$chr} -> {$bin} -> {$line[3]} = $count;
	}
}
close IN;

# count #reads for each gene:

# release count hash
%count = ();

my (%chrgene, $genepos, %geneinfo, @genepos, %genecount);
# count for each gene:
open(IN, $ARGV[1]) || die;
while (<IN>) {
	chomp;
	@line = split(/\t/, $_);
	if ($line[1] ne "Chr") {
		if (/UNKNOWN/) {
			$chr = 0;
		} elsif ($line[1] =~ /(\d+)/) {
			$chr = $1;
		}
		$genepos = $line[2]."\t".$line[3];
		$bin = int($line[2]/BINSIZE);
		$chrgene{$chr} -> {$bin} -> {$line[0]} = $genepos;
		$bin = int($line[3]/BINSIZE);
		$chrgene{$chr} -> {$bin} -> {$line[0]} = $genepos;
		$genecount{$line[0]} = 0;
		$geneinfo{$line[0]} = $chr."\t".$genepos;
	}
}
close IN;

# for each chr:
foreach my $chromosome (sort {$a <=> $b} keys %chrgene) {
	if (exists $chr{$chromosome}) {
		my %bingene = %{$chrgene{$chromosome}}; # bins in this chr
		my %binread = %{$chr{$chromosome}}; # bins in this chr
	#============= start count ================#
	# search all reads in the bin and count:
		foreach my $eachbin (keys %bingene) {
			my %gene = %{$bingene{$eachbin}}; # genes in this bin
			if (exists $binread{$eachbin}) {
			# count #reads for all genes in this chr
				my %pos = %{$binread{$eachbin}}; # all reads in this bin;
				foreach my $pos (keys %pos) {
				# each gene
					foreach my $gene (keys %gene) {
						$genepos = $gene{$gene};
						@genepos = split(/\t/,$genepos);
						if ($pos>=$genepos[0] and $pos<=$genepos[1]) {
							$genecount{$gene} = $genecount{$gene} + $pos{$pos};
						}
					}
				}
			}
		}	# end count
	#============= end count ================#
		# clean the hash:
		delete $chrgene{$chromosome};
		delete $chr{$chromosome};
	} else {
		print STDERR "$chromosome does not exist.\n";
	}
}

# output:
foreach my $genename (keys %genecount) {
	print "$genename\t$geneinfo{$genename}\t$genecount{$genename}\t$geno\n";
}
