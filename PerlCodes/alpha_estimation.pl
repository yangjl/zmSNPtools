#! user/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

my $man = '';
my $help = '';
my $dir = '';
my $output = '';
my $loc_len = '';
## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help,
           man => \$man,
           'dir=s' =>\$dir,
           'output=s' => \$output,
           'length_list=s' => \$loc_len,
          )
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage(1) unless ($dir && $output && $loc_len);


#my $dir = "tajima_D/data/all_dup_QQ_smp_dir_CG";
#my $loc_len = "tajima_D/data/all_dup_QQ_cor_covered_CG";
open(INPUT,"$loc_len");
my %len;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$len{$info[0]}=$info[1];
} 
opendir(DIR,$dir) || die "impossible accedere a $dir";
my @files = grep { -f} map {"$dir/$_"} readdir(DIR);
closedir(DIR);
my $AT_id;
#my $output = "tajima_D/data/all_dup_QQ_cor_CG_alpha";

open(OUTPUT,">$output");
my @total;
my @total_theta;
my $count2=0;
my $max_theta=0;
my $sum_len=0;
my $sum_count=0;
my @sum_theta;
my $sum_theta_count=0;
my $samp_num=0;
my $line_num=0;

for my $file (@files){
	if($file =~ /^$dir\/(\S+)/){
	$AT_id = "$1";
	}
	my $sum_theta1=0;
	my $len = $len{$AT_id};
	$sum_count++;
	open(INPUT,$file);
	my $inside_line_num=0;
        $line_num=0;
	my $flag=0;
	my %seq;
	while(my $line = <INPUT>){
	chomp $line;
	if($line){
	$flag=1;
        }
	my @info = split(/\s+/,$line);
	$samp_num = $#info-2;
	for(my $j=2; $j<=$#info; $j++){
	my $temp=$j-2;
	$seq{$temp}->{$line_num}=$info[$j];
	}

	$line_num++;
	$inside_line_num++;
	}
	if($flag==1){
		for(my $l=$inside_line_num; $l<$len; $l++){
		for(my $h=0; $h<=$samp_num; $h++){
		
		$seq{$h}->{$line_num}="0";
		}
		$line_num++;

		}
	}
	if($flag==1){
	my $temp_seq =print_out($line_num,$samp_num,\%seq); 
	my $seq_alpha= "seq_alpha";
	`perl diff_two_seq.pl $temp_seq $seq_alpha`;		

	my ($alpha,$mean,$var) = get_alpha($seq_alpha);
	print OUTPUT "$AT_id\t$alpha\t$mean\t$var\n";	
	}
	
}

sub get_alpha{
 my $alpha_file =shift(@_);
 open(INPUTa,$alpha_file);
 my @array;
 while(my $linea = <INPUTa>){
 chomp $linea;
 push(@array,$linea);
 }
 my $mean = mean(\@array);
 my @a = map (($_ - $mean)**2, @array);
 my $var= sum(\@a)/scalar @array;
 my $alpha;
 if($var == 0){
 $alpha = 1000000;
 }else{
 if($var > $mean){
 $alpha = ($mean**2)/($var - $mean);
 }else{
  $alpha = ($mean**2)/($var*exp(-$mean/$var));
  }
 }
 return($alpha,$mean,$var);

}

sub sum{
my $arrayref = shift(@_);
my $result;
foreach(@$arrayref){
$result += $_;
}
return $result;
}

sub print_out{
my ($line_num,$samp_num,$hash)= (@_);
my %seq = %{$hash};
my $output_file= "temp_seq_CG";
open(OUTPUT1,">$output_file");
for(my $n=0; $n<=$samp_num; $n++){
print OUTPUT1 "\n$n\t";
	for(my $m=0; $m<$line_num; $m++){
	print OUTPUT1 "$seq{$n}->{$m}";
	}
}
return($output_file);
}

sub mean{
my $array = shift(@_);
my @array = @{$array};
my $sum=0;
my $count1=0;
for(my $j=0; $j<=$#array; $j++){
$sum += $array[$j];
$count1++;
}
my $ave = $sum/$count1;
return($ave);
}

__END__

=head1 NAME

alpha_estimation.pl - estimating alpha parameter of gamma distribution, which describs the distribution of mutation rates among the sites,  from methylation mutation polymorphism. To run this scrip, you need to download diff_two_seq.pl in the same dir, which will be used by alpha_estimation.pl

=head1 SYNOPSIS

Dm_test.pl [options]

 Options:
   -help            brief help message
   -man             full documentation
   -dir             the fold containing the SMP file of each locus
   -output          output file of alpha value
   -length_list     the list of sequence length of mapped methylation state in each locus

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-dir>

The dir is a fold name contaning the SMP file of each locus.

The file format under input dir in Tab delimited text format. For example,

 #chr	position 	sample1	sameple2 sameple3	...
 chr01	1354		0	1	 0		...
 chr01	1547		0	1	 1		...
 ...

=item B<-output>

Output alpha value for each locus(SMP file) to a file.

The output file format. negative_bionomial(NB). For example,

locus_name alpha_value mean_NB variance_NB

=item B<-length_list>

The file format of length_list. For example,

locus_name	the_sequence_length_of_methylation_state

QQS		2000 

Note: locus_name should match the file name under input fold

=back

=head1 DESCRIPTION

B<Alpha_estimation> will estimate alpha value for each locus based on its SMP file

=head1 REFERENCE

A Test for Detecting Selection on DNA Methylation Polymorphisms

=head1 AUTHOR 

Jun Wang @ Wayne State University

=head1 EMAIL

wjeeuc1@gmail.com

=head1 BUGS

none.

=head1 COPYRIGHT 

To request this program, please send an email to wjeeuc1@gmail.com or cfan@wayne.edu 

