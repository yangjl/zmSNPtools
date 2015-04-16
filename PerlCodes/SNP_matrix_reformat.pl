#!/usr/bin/perl -w
# last update:6/23/2014

use strict;
use warnings;

my $file = "$ARGV[0]";
open(IN, $file) || die "Could not open file: $file!\n";

#first line
my $firstline = <IN>;
print "$firstline";

my $chr0 = 0;
my $pos0 = 0;
while (<IN>){
        chomp;
	my @line = split(/,/);
	my $chr = $line[1];
	my $pos = $line[2];
	my $ref = $line[3];
	my $alt = $line[4];

	if($chr == $chr0 && $pos == ($pos0+1)){
		$ref0 = $ref0.$ref;
		$alt0 = $alt0.$alt;
		if($tot > $tot0){
		}	
	}	
	
	
	my @chr = split(/r/, $line[0]);
	print "$chr[1]";
	if ($line[1] < 100000000){
	printf "\tPZE";
	printf '%02s', $chr[1];
	printf '%08s', $line[1];
	}else{
	printf "\tPZE";
        printf '%02s', $chr[1];
        printf '%09s', $line[1];
	}
	
	print "\t$line[1]";

	for my $i (2..$#line){
		if($line[$i] =~ /--/){
		print "\tN";
		}else{
		my @geno = split(/\//, $line[$i]);
		$B73 = $geno[0];
		print "\t$geno[1]";
		}
	}
	print "\t$B73";#B73

        print "\n";
}
close IN;

