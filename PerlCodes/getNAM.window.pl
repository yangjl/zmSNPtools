#!/usr/bin/perl -w
#Jinliang Yang
#12/3/2009
#This script is used to extract the slide window of 10K. 

use strict;
use warnings;

#the output file is Hp301toFGS
my $out= $ARGV[1];

open (OUTFILE, ">$out") or die "Can not open the output file '$out'$!\n";

#open Bowtie aligned file:
# load all reads info into hash
open(IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	@line = split(/\t/, $_);
	if (/chr(\d+)/) {
		$chr = $1;
	}
	if ($chr =10){
		if ($line[3] >= 83500000 and $bline[2] <= 87500000){	
			print OUTFILE "@line\n";
		}
	}
}
close IN;
