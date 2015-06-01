#!/usr/bin/perl -w
#Jinliang Yang
use strict;
use warnings;

my $input = $ARGV[0];
my $output = $ARGV[1];

my (@bline, $winstart, $winend, %count);
my $slidewin =10000;

open(IN, $input) || die;


while(<IN>){
	chomp;
	@bline =split(/,/, $_);
	$winstart = 83500000;
	$winend =83500000+$slidewin;
	if($bline[1]>=$winstart and $bline[1]<=$winend){
		$count{$winstart}++;
		}
	$winstart = $winstart + $slidewin;
	$winend = $winend + $slidewin;
}
open(OUTF, $output) || die;
print OUTF "$winstart,$count";
close(IN);
close(OUTF);