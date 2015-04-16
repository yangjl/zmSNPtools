#!/usr/bin/perl -w
# usage: perl INDEL_merge.pl [input] > [output]
# last update: 6.23.2014

use strict;
use warnings;

my ($count, $chr0, $pos0) = (0, "1", 0);
my @linearray;
my $file = "$ARGV[0]";

open(IN, $file) || die "I Could not open the INDEL matrix: $file!\n";

while (<IN>){
    chomp;
    if($_ !~ /^#/ && $_ =~ /\*/) {
    	$count++;
    	
        my @line = split(/,/);
        my $chr = $line[1];
        my $pos = $line[2];
        
        if($count == 1){
            push @linearray, $_;
            $chr0 = $chr;
            $pos0 = $pos;
        } elsif($count > 1 && $chr == $chr0 && $pos == $pos0+1){
            push @linearray, $_;
            $chr0 = $chr;
            $pos0 = $pos;
        } else{
            $chr0 = $chr;
            $pos0 = $pos;
            
            #### print the indel chunk, use the first SNPs' id, chr and pos
            my @tem = split(/,/, $linearray[0]);
            print "$tem[0],$tem[1],$tem[2]";
            
            ### join the indel bps and get line of the most informative base!
            my $totcount = 0;
            my ($temref, $temalt, $wmax, $temcount);
            for my $x (0..$#linearray){
                my @temline = split(/,/, $linearray[$x]);
                $temref .= $temline[3];
                $temalt .= $temline[4];
                
		$temcount =0;
                for my $y (5..$#temline){
                    $temcount += $temline[$y];
                }
		#print "\t", $temcount;
                if($temcount > $totcount){
                	$totcount = $temcount;
                	$wmax = $x;
                }
                
            }
            
            print ",$temref,$temalt";
            
            ### print the most informative base!
            my @temarray = split(/,/, $linearray[$wmax]);
            for my $i (5..$#temarray) {
                print ",", $temarray[$i];
            }
            print "\n";
            
            @linearray = ();
            push @linearray, $_;
        }
	}
}



