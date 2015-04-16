#!/usr/bin/perl -w
# Jinliang Yang
# Version: Diallel version1 for GenSel and PLINK 
# Usage: ./myimpute_diallel.pl dsnp=dsnp_4impute[dsnp set] diallel=diallel_ped.txt mode=1 > output
# update: 4/17/2012
# 

use strict;
use warnings;
use Getopt::Long;

my $diallel = "diallel_4impute.txt";
my $dsnp = "fdsnp_2M_012612.dsnp";
my $mode = 1; #1 for GenSel; #2 for PLINK
my $header = 1; # 1 with header; #2 without header

GetOptions (
            'dsnp=s' => \$dsnp, #density SNP 
            'diallel=s' => \$diallel,
            'mode=i' => \$mode,
            'header=i' => \$header);


#### Read in the diallel parents
open(FILE, $diallel) || die "I could not open genotype of diallel file: $diallel!\n";
my $col_names =<FILE>;
my @parents;
while (<FILE>){
    chomp;
    my @line = split(/\t/, $_);
    my $cross = $line[2];
    push @parents, $cross;
}
close FILE;

###### Read in the high density SNP file:
open(IN, $dsnp) || die "I Could not open the parental dSNP file: $dsnp!\n";
my $first_line = <IN>;
chomp $first_line;
my @names =split(/\t/, $first_line);

my @snp2M;
while (<IN>) {
    chomp;
    my @line = split(/\t/, $_);
    
    #### create a "array of hash" for each snp  
    my $snptype;
                
        for (my $i=0; $i<=$#line; $i++){
            $snptype -> {$names[$i]} = $line[$i];
            }# End of the for loop
        push @snp2M, $snptype;
}
close IN;

###### coding system ###############
my ($BB, $AA, $AB, $BN, $AN, $NN);
if($mode==1){
        $BB=-10;
        $AA=10;
        $AB=0;
        $BN=-5;
        $AN=5;
        $NN=0;
}else{
        $BB="B B";
        $AA="A A";
        $AB="A B";
        $BN="N N";
        $AN="N N";
        $NN="N N";
}


#### Print the header line
if($header==1){
    print "ID";
    for my $snpid (0..$#snp2M){
        print "\t$snp2M[$snpid]{'rs'}";
    }
    print "\n";
}

#### Impute the parental genotypes
for my $f1 (0..$#parents){
    my @tem = split(/x/, $parents[$f1]);
    my $p1 = $tem[0];
    my $p2 = $tem[1];
    
    if($mode==1){
        print "$parents[$f1]";
    }else{
        print "$parents[$f1]\t1\t$p1\t$p2\t1";
    }

    my $count=0;
    for my $snpid (0..$#snp2M){
        my $snpp1 = $snp2M[$snpid]{$p1};
        my $snpp2 = $snp2M[$snpid]{$p2};
        my $snpb73 = $snp2M[$snpid]{'B73'};
        my $geno;
        $count++;

        if($snpp1 eq $snpb73){
            if($snpp2 ne "N" && $snpp2 ne $snpb73){
                $geno=$AB;
            } elsif ($snpp2 eq $snpb73){
                $geno=$BB;
            } elsif ($snpp2 eq "N"){
                $geno=$BN;
            }
        }elsif($snpp1 eq "N"){
            if($snpp2 ne "N" && $snpp2 ne $snpb73){
                $geno=$AN;
            }elsif($snpp2 eq "N"){
                $geno=$NN;
            }elsif($snpp2 eq $snpb73){
                $geno=$BN;
            } 
        }elsif($snpp1 ne "N" && $snpp1 ne $snpb73){
            if($snpp2 ne "N" && $snpp2 ne $snpb73){
                $geno=$AA;
            }elsif($snpp2 eq $snpb73){
                $geno=$AB;
            }elsif($snpp2 eq "N"){
                $geno=$AN;
            }
        }else{$geno=$NN;}
        print "\t$geno";
     }
     print "\n";
     print STDERR "###---Finished imputing $count SNP for Diallel $parents[$f1]!---###\n";
 }




