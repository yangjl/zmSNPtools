# !usr/bin/perl -w
# Jinliang Yang 11-29-2009

use strict;
use warnings;
use Tie::File;

my $filename = 'FGS.part';
tie my @FGS, 'Tie::File', $filename or die "Cannot tie '$filename' $!\n";

my $nam = 'Hp301.part';
tie my @NAM, 'Tie::File', $nam or die "Cannot tie '$nam' $!\n";

#the output file is Hp301toFGS
my $out= 'Hp301toFGS.test';
open(OUTFILE, ">$out") or die "Can not open the output file '$out'$!\n";

my $i =0;

	#FGS.part
	#1.chr0    GRMZM2G056896   GRMZM2G056896_T01       GRMZM2G056896_P01       yes     216025  216655  -1      cdna|estgood
	foreach(@FGS){
		my($f1, $f2, $f3, $f4, $f5, $f6, $f7, $f8, $f9, $f10) = split(/\t/, $_);
		while(<@NAM>){
		chomp $_;	
		#Hp301.part
		#SRR026759.4.1 HWI-EAS210_Buckler_71221_7_1_111_118.1 length=36  -       3.chr1    4.300049167       CCATTTCCTAGCCGTTG
		#CTCCATTTCTGCCGACCGC     +*!+()+++!+++!++"+&+++++++++++++++++    0	
		my($n1, $n2, $n3, $n4, $n5, $n6, $n7) = split(/\t/, $_);	
		if ($f1 =~ /$n3/){
			if (($f6 >= $n4 and $f7 <= $n4) || ($f7 >= $n4 and $f6 <= $n4)) { 
				$i++;
				next;
			}
		}
	}
			print OUTFILE "$f1\t$f2\t$f6\t$f7\t$i\n";

}
close (OUTFILE);
exit;



