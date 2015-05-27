#!/usr/local/bin/perl -w
use Getopt::Long;
use Pod::Usage;

my $man = '';
my $help = '';
my $input = '';
my $output = '';
my $windowsize = '';
my $alpha = '1';
## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help,
			man => \$man,
			'input=s' => \$input,
			'output=s' => \$output,
			'length=i' => \$windowsize,
			'alpha=f' =>\$alpha,
			)
or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage(1) unless ($input && $output && $windowsize);



open (INPUT,"$input") or die "can't open input FILE: $!";
while(<INPUT>)
{
	chop;
	/^#/ and next;
	$_ =~ s/\s*$//g;
	@line = split (/\t/,$_);
	$chr=shift @line;
#	$position=shift @line;
	shift @line;
	$samplenumber=$#line+1;
	print "$samplenumber\n";
	if($samplenumber < 2)
	{
		print "The amount of sample is less than 2, please check your data.\n";
		exit;
	}
	%genotype=();
	foreach (@line)
	{	#if($_){
		if(!exists $genotype{uc($_)})
		{
			$genotype{uc($_)}=1;
		}
		else
		{
			$genotype{uc($_)}++;
		}
		#}
	}
	my $key_num = keys %genotype;
	next if(keys %genotype == 1);
	$minnumber=$samplenumber;
	#print "$minnumber\t";
	foreach (keys %genotype)
	{
		$minnumber = $genotype{$_} if($genotype{$_} < $minnumber);
	#	print "min:$minnumber\n";
	}
#	$potentialminbinnumber=int (($position-1)/($windowsize));
#	$potentialminbinnumber=0 if $potentialminbinnumber < 0;
#	$potentialmaxbinnumber=int (($position-1)/($windowsize-0));
        #print "$potentialminbinnumber\t$potentialmaxbinnumber\n";
#	for ($binnumber=$potentialminbinnumber;$binnumber<=$potentialmaxbinnumber;$binnumber++)
#	{	
		if (!exists $snp{$chr})
		{
			$snp{$chr}=1;
		}
		else
		{
			$snp{$chr}++;
		}
		if (!exists $pisnp{$chr})
		{
			$pisnp{$chr}{$minnumber}=1;
		}
		else
		{
			$pisnp{$chr}{$minnumber}++;
		}
#	}
	
}
close (INPUT);
$a1=0;
$a2=0;
for($i=1;$i<$samplenumber;$i++)
{
	$a1 += 1/$i;
	$a2 += 1/($i*$i);
}
$a3 = ($a1*$a1 - $a2)/2;
$b1=($samplenumber+1)/(3*($samplenumber-1));
$b2=2*($samplenumber*$samplenumber+$samplenumber+3)/(9*$samplenumber*($samplenumber-1));
$c1 = 2*$a1 - (3*$a3/$a1);
$c3 = $b1 - 1/$a1;
$c4 = $b2 - ($samplenumber+2)/($a1*$samplenumber)+$a2/($a1*$a1);
$e1 = $c3/$a1;
$e2 = $c4/($a1*$a1+$a2);
print "a1:$a1\ta2:$a2\ta3:$a3\tb1:$b1\tb2:$b2\tc1:$c1\tc3:$c3\tc4:$c4\te1:$e1\te2:$e2\n";
open (OUTPUT,">$output") or die "can't open OUT-FILE: $!";

	print OUTPUT "#chr\tstart\tend\tDm\tsegregation_site\ttheta_pi\ttheta_s\n";
foreach $chr (sort keys %snp)
{
#	foreach $bin (sort {$a <=> $b} keys %{$snp{$chr}})
#	{	
		$binstart=1;
		$binend=$windowsize;
		$s=$snp{$chr};
		$theta=$s/($a1*$windowsize);
		$sumx=0;
		foreach (keys %{$pisnp{$chr}})
		{
			$sumx +=$_*($samplenumber-$_)*$pisnp{$chr}{$_};
		}
		
		$pi=$sumx*2/($samplenumber*($samplenumber-1)*$windowsize);
#		$td= ($sumx*2/($samplenumber*($samplenumber-1))-$s/$a1)/sqrt($e1*$s+$e2*$s*($s-1));
		$theta_pi = $pi * exp(2*($alpha +1)*$pi/($alpha));
		$theta_s = $theta * exp($c1*($alpha+1)*$theta/$alpha);
		$wd = ($theta_pi - $theta_s)/sqrt($a1*$e1*$theta_s/$windowsize+$e2*$a1*$theta_s*($a1*$theta_s-1/$windowsize));
			print OUTPUT "$chr\t$binstart\t$binend\t$wd\t$s\t$theta_pi\t$theta_s\n";

#	}
}
close (OUTPUT);
__END__

=head1 NAME

Dm_test.pl - Creating Dm test value for methylation mutation polymorphism

=head1 SYNOPSIS

Dm_test.pl [options]

 Options:
   -help            brief help message
   -man             full documentation
   -input           input file
   -output          output file
   -length          the sequence length of methylation state
   -alpha           the alpha parameter of gamma function for the distribution of mutation rate, default is 1

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-input>

Input data from a file in Tab delimited text format. For example,

 #chr	position 	sample1	sameple2 sameple3	...
 chr01	1354		0	1	 0		...
 chr01	1547		0	1	 1		...
 ...

=item B<-output>

Output data to a file.

=item B<-length>

Set the sequence length of methylation state

=item B<-alpha> 

the alpha value will determine the gamma distribution which describs the frequency distribution of mutation rate. 
The smaller the alpha value, the more variable the mutation rates among the sites

=back

=head1 DESCRIPTION

B<Dm_test> will calculate Dm test value.

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

