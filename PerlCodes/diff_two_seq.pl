#! user/bin/perl -w
use strict;
my $seq_file  =$ARGV[0];
my $output = $ARGV[1];

my %total_dist;
open(INPUT,"$seq_file");
<INPUT>;
my %seq;
my @seq;
my $seq_count=0;
my $seq_len=0;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
my @seq_sub = split(//,$info[1]);
my $seq="";
$seq = $info[1];
$seq_len = length($seq);

$seq{$seq_count}=$seq;
$seq_count++;
}
my $num = (keys %seq )-1;
my %dist;
for(my $i=0; $i <=($num-1); $i++){
my $seq1= $seq{$i};
	for(my $k=$i+1;$k<=$num; $k++){
	my $seq2 = $seq{$k};
	my $dist = seq_dist($seq1,$seq2);
	my $key = "$i"."_"."$k";
	$total_dist{$key}=$dist;
	$dist{$key} = $dist;
	}
}
my %already;
my @already;
my $key = small_dist_key(\%dist);
my @index = split(/_/,$key);

for(my $j=0; $j<=$#index; $j++){
push(@already,$index[$j]);
$already{$index[$j]}=1;
}

my $rest =$num -2;
for(my $h=0; $h<=$rest; $h++){ 
my $sum_dist=0;
my %join_dist;
for(my $n=0; $n<=$num; $n++){
if(!(exists $already{$n})){
for(my $h=0; $h<=$#already;$h++){
my $temp_key= "$already[$h]"."_"."$n";
my $temp_key1 = "$n"."_"."$already[$h]";
my $dist;

if(exists $total_dist{$temp_key}){
$dist = $total_dist{$temp_key};
}elsif(exists $total_dist{$temp_key1}){
$dist = $total_dist{$temp_key1};
}else{
my $already_seq_id = $already[$h];
my $already_seq = $seq{$already_seq_id};

my $compare_seq = $seq{$n};
$dist = seq_dist($already_seq,$compare_seq);
my $temp_key2 = "$already_seq_id"."_"."$n";
$total_dist{$temp_key2}=$dist;
}
$sum_dist += $dist;
}
$sum_dist = $sum_dist/($#already + 1);
$join_dist{$n} = $sum_dist;
}
}

my $small_dist_key = small_dist_key(\%join_dist);
push(@already,$small_dist_key);
$already{$small_dist_key}=1;
}


my %tree_seq;
for (my $m=0; $m<=$num; $m++){
my $tmp_seq=$seq{$already[$m]};
my @tmp_seq = split(//,$tmp_seq);
for(my $m1=0; $m1<=$#tmp_seq; $m1++){
$tree_seq{$m}->{$m1}=$tmp_seq[$m1];

}
}

my $change_num=0;
my $sub_0=0;
my $sub_1=0;
my $sub_2=0;
my $sub_3=0;
my $sub_4=0;
my $sub_5=0;
my $sub_6=0;
for (my $h=0; $h<=$seq_len-1; $h++){
my %genotype;
my %change_time;
$change_time{count}=0;
$genotype{$tree_seq{0}->{$h}}++;
my $pre_status=$tree_seq{0}->{$h};
for(my $h1=1; $h1<=$num; $h1++){
$change_time{$change_time{count}}++;

if($tree_seq{$h1}->{$h} != $pre_status){
$change_time{count}++;
$pre_status = $tree_seq{$h1}->{$h};
}
$genotype{$tree_seq{$h1}->{$h}}++;

}
$change_time{$change_time{count}}++;

my $key_num = keys %genotype;
if($key_num==1){
$sub_0++;
}else{
if(($genotype{0}==1)||($genotype{1}==1)){
$sub_1++;
}elsif($change_time{count}==1){
$sub_1++;
}elsif(($genotype{0}==2) || ($genotype{1}==2)){
$sub_2++;
}elsif($change_time{count}==2){
$sub_2++;
}elsif(($change_time{count}==3)&&(($change_time{1}==1) || ($change_time{2}==1))){
$sub_2++;
}elsif(($genotype{0}==3) || ($genotype{1}==3)){
$sub_3++;
}elsif($change_time{count}==3){
$sub_3++;
}elsif(($change_time{count}==4) && (($change_time{1}==1) || ($change_time{2}==1) || ($change_time{3}==1))){
$sub_3++;
}elsif(($change_time{count}==5)&&((($change_time{1}==1)&&($change_time{3}==1))||(($change_time{2}==1)&&($change_time{4}==1)))){
$sub_3++;
}elsif(($genotype{0}==4)||($genotype{1}==4)){
$sub_4++;
}elsif($change_time{count}==4){
$sub_4++;
}elsif(($change_time{count}==5)&&(($change_time{1}==1)||($change_time{2}==1)||($change_time{3}==1)||($change_time{4}==1))){
$sub_4++;
}elsif(($change_time{count}==6)&&((($change_time{1}==1)&&($change_time{3}==1))||(($change_time{2}==1)&&($change_time{4}==1))||(($change_time{3}==1)&&($change_time{5}==1)))){
$sub_4++;
}elsif(($change_time{count}==7)&&((($change_time{1}==1)&&($change_time{3}==1)&&($change_time{5}==1)) ||(($change_time{2}==1)&&($change_time{4}==1)&&($change_time{6}==1)))){
$sub_4++ 

}elsif(($genotype{0}==5)||($genotype{1}==5)){
$sub_5++;
}elsif($change_time{count}==5){
$sub_5++;
}elsif(($change_time{count}==6)&&(($change_time{1}==1)||($change_time{2}==1)||($change_time{3}==1)||($change_time{4}==1)||($change_time{5}==1))){
$sub_5++;
}elsif(($change_time{count}==7)&&((($change_time{1}==1)&&($change_time{3}==1))||(($change_time{2}==1)&&($change_time{4}==1))||(($change_time{3}==1)&&($change_time{5}==1))||(($change_time{4}==1)&&($change_time{6}==1)))){
$sub_5++;
}elsif(($change_time{count}==8)&&((($change_time{1}==1)&&($change_time{3}==1)&&($change_time{5}==1))||(($change_time{2}==1)&&($change_time{4}==1)&&($change_time{6}==1))|| (($change_time{3}==1)&&($change_time{5}==1)&&($change_time{7}==1)))){
$sub_5++;
}elsif(($change_time{count}==9)&&((($change_time{1}==1)&&($change_time{3}==1)&&($change_time{5}==1)&&($change_time{7}==1))||(($change_time{2}==1)&&($change_time{4}==1)&&($change_time{6}==1)&&($change_time{8}==1)))){
$sub_5++;
}else{
$sub_6++;
}
}
}

open(OUTPUT,">$output");
for(my $t=0; $t<($sub_0);$t++){
print OUTPUT "0\n";
}
for(my $t=0; $t<$sub_1; $t++){
print OUTPUT "1\n";
}
for(my $t=0; $t<$sub_2; $t++){
print OUTPUT "2\n";
}
for(my $t=0; $t<$sub_3; $t++){
print OUTPUT "3\n";
}
for(my $t=0; $t<$sub_4; $t++){
print OUTPUT "4\n";
}
for(my $t=0; $t<$sub_5; $t++){
print OUTPUT "5\n";
}
for(my $t=0; $t<$sub_6; $t++){
print OUTPUT "6\n";
}
sub small_dist_key{
my $hash = shift(@_);
my %hash = %{$hash};
my $small_key="";
for my $key (sort { $hash{$a} <=> $hash{$b}} keys %hash){

$small_key = $key;
goto line1;
}
line1:return($small_key);

}



sub seq_dist{
my ($seq1,$seq2) = (@_);
my @info1 = split(//,$seq1);
my @info2 = split(//,$seq2);
my $len = $#info1+1;
my $diff=0;

	for(my $n=0; $n<=$#info1; $n++){
		if($info1[$n] ne $info2[$n]){
			$diff++;
		}
	}
#	my $dist = -3*log(1-4*$diff/($len*3))/4;
#	my $dist = -log(1-2*$diff/($len))/2;
	#print "$dist\n";
	#exit;
	my $dist= $diff;
return($dist);
}






