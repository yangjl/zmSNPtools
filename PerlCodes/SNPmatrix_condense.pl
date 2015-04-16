#!/usr/bin/perl -w
# usage: perl INDEL_merge.pl [input] > [output]
# last update: 6.23.2014
# Note input file could be ordered by chr and pos!

use strict;
use warnings;

my ($count, $chr0, $pos0) = (0, "1", 0);
my @linearray;
my $file = "$ARGV[0]";

open(IN, $file) || die "I Could not open the INDEL matrix: $file!\n";

while (<IN>){
    chomp;
    #### output SNP loci
    if($_ !~ /^#/ && $_ !~ /\*/){
    	my @line = split(/\t/);
    	my $chr = $line[0];
    	my $pos = $line[1];
    	$chr =~ s/chr//g;
    	print "$chr\_$pos\t$chr\t$pos\t$line[2]";
    	for my $i (3..$#line){
    		if($line[$i] eq "--"){
    			print "\tN";
    		}else{
    			my @cell = split(/\//, $line[$i]);
    			print "\t$cell[1]";
    		}	
    	}
    	print "\n";
    }
    #### output Indel loci
    elsif($_ !~ /^#/ && $_ =~ /\*/) {
    	$count++;
    	
        my @line = split(/\t/);
        my $chr = $line[0];
        my $pos = $line[1];
        
        if($count == 1){
            push @linearray, $_;
            $chr0 = $chr;
            $pos0 = $pos;
        } elsif($count > 1 && $chr eq $chr0 && $pos == $pos0+1){
            push @linearray, $_;
            $chr0 = $chr;
            $pos0 = $pos;
        } elsif($count > 1 && $chr eq $chr0 && $pos < $pos0+1) {
        	die "The file should be ordered by chr and pos: $file!\n"
        } else {
            $chr0 = $chr;
            $pos0 = $pos;
            
            #### print the indel chunk, use the first SNPs' id, chr and pos
            my @tem = split(/\t/, $linearray[0]);
            $tem[0] =~ s/chr//g; # subsititution
            print "$tem[0]\_$tem[1]\t$tem[0]\t$tem[1]";
            
            ### select the most informative base!
            my $totmiss = 99999;
            my ($temref, $wmax);
            for my $x (0..$#linearray){
                my @temline = split(/\t/, $linearray[$x]);
                #$temref .= $temline[3];
                #$temalt .= $temline[4];
                my $temmiss = 0;
                for my $y (3..$#temline){
                	if($temline[$y] =~ "--"){
                		$temmiss += 1;
                	}
                    
                }
                if($temmiss < $totmiss){
                	$totmiss = $temmiss;
                	$wmax = $x;
                }    
            }
            ### get the ref allele info
            my @maxline = split(/\t/, $linearray[$wmax]);
            for my $i (3..$#maxline){
            	if($maxline[$i] !~ "--"){
            		my @alleles = split(/\//, $maxline[$i]);
                	if($alleles[0] =~ /\*/){
                		$temref = "-";
                	} elsif($alleles[1] =~ /\*/){
                		$temref = "+";
                	}
                } 
            }
            print "\t$temref";
            
            ###### print the indels
            for my $i (3..$#maxline){
            	if($maxline[$i] =~ "--"){
            		print "\tN";
            	}else{
            		my @alleles = split(/\//, $maxline[$i]);
            		if($alleles[0] eq $alleles[1]){
            			print "\t$temref";
            		}elsif($alleles[0] =~ /\*/){
            			print "\t+";
            		}elsif($alleles[1] =~ /\*/){
            			print "\t-";
            		}
            	}
            }
            print "\n";
            
            @linearray = ();
            push @linearray, $_;
        }
    }
}


