#!/usr/bin/perl -w
# Jinliang Yang
# Version: version3(beta) for GenSel 
# Usage: ./myimpute_beta_4GenSel.pl dsnp=dsnp_4impute[dsnp set] tsnp=tsnp_4impute.txt[tsnp set] > output

#####################################################################
# last update:03/07/2011
# Modification: 
# 1. Generate two versions: one is for PLINK use, the other is for GenSel use
# 2. Can handle the SNP in one file
#
#####################################################################

use strict;
use warnings;
use Getopt::Long;

my $dsnp = "dsnp_4impute";
my $tsnp = "tsnp_4impute.txt";

my $header = 1; #print the header; 2 no header
my $mode = 1; #1 for GenSel; 2 for PLINK
my $physical = 1; #1 using physical position, 2 using genetic position
my $start = 6; # start of the RIL #
my $end = 4897; # End of the RIL #
my $prob = 0; # probability used as the threshold for SNP projection


GetOptions (	
				'dsnp=s' => \$dsnp, #density SNP 
				'tsnp=s' => \$tsnp, #tag SNP
				
				'header=i' => \$header, #1 header, #2 no header
				'mode=i' => \$mode, #1 for GenSel #2 for PLINK.
                'physical=i' => \$physical, #1 for physical position[pos], 2 for genetics pos [genetpos]
                'start=i' => \$start,
                'end=i' => \$end,
                'prob=i' => \$prob);



######-----------------------------------------------
#
# Part1: read in the parental snp data:
#
#####------------------------------------------------
open(IN, $dsnp) || die "I Could not open the parental dSNP file: $dsnp!\n";
my $first_line = <IN>;
chomp $first_line;
my @names =split(/\t/, $first_line);
my @snp;

while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	
	#### create a "array of hash" for each snp  
	my $snptype;
	
	 	for (my $i=0; $i<=$#line; $i++){
		$snptype -> {$names[$i]} = $line[$i];	
		}# End of the for loop
		push @snp, $snptype;
}
close IN;

#End of printing the first line:
if ($header==1){
# Print The SNP ID:
		print "ID";
		for my $j (0..$#snp){
			my $temp = $j+1;
			print "\t$snp[$j]{'chr'}" ."_". "$snp[$j]{'pos'}";
			}
		print "\n"; 
}


######------------------------------------------------
#
# Part2: read in the 5000 NAM RIL 1000 markers data:
#
######------------------------------------------------

open(MAP, $tsnp) || die "I could not open genetic map tSNP file: $tsnp!\n";
my @loci; ## array of array
my $colcount;
while (<MAP>) {
	chomp;
	my @line = split(/\t/, $_);
	$colcount= $#line;
	
	my $locus;
	#### create a "array of array" for each locus:
	
	    
		for (my $i=0; $i<=$#line; $i++){	
			
		#for (my $i=0; $i<=100; $i++){
		$locus -> [$i] = $line[$i];	
		}
		push @loci, $locus; 
}
close MAP; # End of this part.

#test code:
#  for my $i (0..$#loci){
#	   print "Row $i is ---------->";
#	    for my $j (0..$#{$loci[$i]}){
#		   print "$loci[$i][$j]\t";
#	   }
#	   print "\n";
#    }

######--------------------------------------------------------------#
#
# Part3: Calculate the NAM RIL break points and project the genotype:
#
######--------------------------------------------------------------#
	

#the dimension should be changed here $start and $end
for (my $col=$start; $col <= $end; $col++){
#for (my $col=4890; $col <= $colcount; $col++){  # first layer of col loop, loop the NAM RIL individually.
	
	#1: pull out the NAM RIL Genotype one by one
	my @rilgeno;
	
	push @rilgeno, [@{$loci[0]}[0..4,$col]];
	for (my $x=1; $x<=$#loci; $x++){
		if ($loci[$x][$col] != -9) {
		push @rilgeno, [@{$loci[$x]}[0..4,$col]];
		}
	}#End of the for loop

#test code:
#  for my $i (0..$#rilgeno){
#	   print "Row $i is ---------->";
#	    for my $j (0..$#{$rilgeno[$i]}){
#		   print "$rilgeno[$i][$j]\t";
#	   }
#	   print "\n";
#    }

	#2. get the ril interval info
	
	my @interval=();
	for my $x (0..$#rilgeno){
		for my $y (0..$#{$rilgeno[$x]}){
		$interval[$x][$y] = $rilgeno[$x][$y];
		}
	}
	$interval[0][6] = "NPos";
	$interval[0][7] = "NGeno";
	$interval[0][8] = "NChr"; # Define another chomosome #
	
	for my $m (1..($#rilgeno-1)){
		$interval[$m][6]=$rilgeno[$m+1][4]; # physical position
		$interval[$m][7]=$rilgeno[$m+1][5]; # genotype 2
		$interval[$m][8]=$rilgeno[$m+1][2]; # Chr # 2
	}
		$interval[$#rilgeno][6]= -9;
		$interval[$#rilgeno][7]= -9;
		$interval[$#rilgeno][8]= -9;
		
 
#test the @interval: 
# 		for my $i (0..$#interval){
#	   print "Row $i is ---------->";
#	    for my $j (0..$#{$interval[$i]}){
#		   print "$interval[$i][$j]\t";
#	   }
#	   print "\n";	
#    }# End of the test

		
    #3:	get the break points of the interval:
		my @break=();	
		for my $row (0 .. $#interval){		
			if (($interval[$row][5] ne $interval[$row][7]) or ($interval[$row][2] ne $interval[$row][8])){
				push @break, [ @{$interval[$row]}[0..8] ];
				} 
				#end of if statment
		}# End of the inner for loop
 	
#   test the @interval: 
# 		for my $i (0..$#break){
#	   print "Row $i is ---------->";
#	    for my $j (0..$#{$break[$i]}){
#		   print "$break[$i][$j]\t";
#	   }
#	   print "\n";	
#    }# End of the test

#print "\n";

		#4: Dereference the break file:
		my $rilname = $break[0][5];
		my @namname = split(/E/, $rilname);
					
			#4.0: Project the RIL name:
			if ($mode==1){
			print "$rilname";
			}else{
			print "$namname[0]\t$namname[1]\t$namname[0]\tB73";
			}
			
			
			
		my $snpchr;
		my $snppos;
		############------------- projection code:-------##############
		my ($B73, $nonB73, $missing);
		if($mode==1){
			$B73=-10;
			$nonB73=10;
			$missing=0;
		}else{
			$B73="B B";
			$nonB73="A A";
			$missing="N N";
		}
		############------------------------------------##############
		
		for my $i (0..$#snp){ # SNP for loop
				my $tem;
				$snpchr = $snp[$i]{'chr'};
				
				if ($physical==1){
				$snppos = $snp[$i]{'pos'};
				} else{
				$snppos = $snp[$i]{'genetpos'};
				}
				
				
				#case one: chr 1:
					if ($snpchr == 1 && $snppos <= $break[1][4]) {
						#---------------------------------------------------------------------------
						if($break[1][5] == 0){
						$tem= $B73; ###$snp[$i]{B73};
						}else{
								####-----------------------------------------------------------------------------
									if($snp[$i]{$namname[0]} eq $snp[$i]{B73} and $snp[$i]{$namname[0]} ne "N" ){
									$tem= $B73;
									} elsif ($snp[$i]{$namname[0]} eq "N"){
									$tem= $missing;
									} else{
									$tem= $nonB73;
									}# end of block
								####-----------------------------------------------------------------------------
						}
					}
					
				
				##############	
				#OTHER CASES #
				##############
					else {	# start of else 				
							for my $t (1..$#break){ # $t for loop
							my $t1=$t+1;
								
					#case two: end of each chromosome.
							if ($snpchr == $break[$t][2] && $snpchr != $break[$t][8] && $snppos >= $break[$t][4]) {
								#---------------------------------------------------------------------------
								if($break[$t][5] == 0){
								$tem= $B73; ###$snp[$i]{B73};
								}else{
								####-----------------------------------------------------------------------------
									if($snp[$i]{$namname[0]} eq $snp[$i]{B73} and $snp[$i]{$namname[0]} ne "N" ){
									$tem= $B73;
									} elsif ($snp[$i]{$namname[0]} eq "N"){
									$tem= $missing;
									} else{
									$tem= $nonB73;
									}# end of block
								####-----------------------------------------------------------------------------
								}
							}
					#case three: start of each chromosome.
							elsif ($snpchr != $break[$t][2] && $snpchr == $break[$t][8] && $snppos <= $break[$t][6]) {
								#---------------------------------------------------------------------------
								if($break[$t][7] == 0){
								$tem= $B73; ###$snp[$i]{B73};
								}else{
									####-----------------------------------------------------------------------------
									if($snp[$i]{$namname[0]} eq $snp[$i]{B73} and $snp[$i]{$namname[0]} ne "N" ){
									$tem= $B73;
									} elsif ($snp[$i]{$namname[0]} eq "N"){
									$tem= $missing;
									} else{
									$tem= $nonB73;
									}# end of block
								####-----------------------------------------------------------------------------
								}
							}
					#case four: non recombinant intervals.
							elsif ($snpchr == $break[$t][8] && $snpchr== $break[$t1][2] && 
									$snppos >= $break[$t][6] && $snppos <= $break[$t1][4]){
								#---------------------------------------------------------------------------
								if($break[$t][7] == 0){
								$tem= $B73; ###$snp[$i]{B73};
								}else{
								####-----------------------------------------------------------------------------
									if($snp[$i]{$namname[0]} eq $snp[$i]{B73} and $snp[$i]{$namname[0]} ne "N" ){
									$tem= $B73;
									} elsif ($snp[$i]{$namname[0]} eq "N"){
									$tem= $missing;
									} else{
									$tem= $nonB73;
									}# end of block
								####-----------------------------------------------------------------------------
								}	
							}
					#case five: Recombinant intervals.
							elsif ($snpchr == $break[$t][2] && $snpchr == $break[$t][8] && 
									$snppos > $break[$t][4] && $snppos < $break[$t][6]){
								#---------------------------------------------------------------------------
								#!!!!! Recombination region:
								# start of else block
								if($snp[$i]{$namname[0]} eq $snp[$i]{B73} and $snp[$i]{$namname[0]} ne "N" ){
									$tem=$B73;
								} elsif ($snp[$i]{$namname[0]} eq "N"){
									$tem= $missing;
								} else { #else0
								####################################################################
								if ($physical==1){
								
								my $dist = 20 * ($snppos - $break[$t][4])/($break[$t][6] - $break[$t][4]);
									if ($break[$t][5] == 0){
									$tem = -10 + $dist;
									} else{
									$tem = 10 -$dist;
									} 
									
									if($mode==1){
									$tem = int($tem);
									}else{ #else3
										if($tem>0){$tem ="A A";}
											elsif($tem==0){$tem = "N N";}
											else{$tem="B B";}
										}#else3
									
									} else { #else2
									
							my $prb1= (1-($snppos-$break[$t][4])/100) * ($break[$t][6]-$snppos)/($break[$t][6]-$break[$t][4]);
							my $prb2=($snppos-$break[$t][4]) * (1- ($break[$t][6]-$snppos)/100)/($break[$t][6]-$break[$t][4]);
							
							if ($prb1 > $prb2 and $prb1 >= $prob){
										if($break[$t][5] == 0){$tem = $B73;}else{$tem=$nonB73;}								
							}elsif ($prb1 < $prb2 and $prb2 >= $prob){
										if($break[$t][5] == 2){$tem=$nonB73;}else{$tem=$B73;}
							}else{$tem=$missing;} 


								} #else2
								} #inside else0	
								####################################################################
							} # of case five else
					}# end of inside for loop
				##############	
				}# END of OTHER CASES #
				##############
									
			
			print "\t$tem";	
			}# End of SNP for loop		
			
			print "\n";

}# End of the Column loop

#############################################################################################################

#############################################################################################################
