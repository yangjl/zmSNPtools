# !usr/bin/perl -w
# Get the DNA sequence data
print "please type file name :~/Centra-SeqDB/Psedudo/V1.20090320/maize_pseudo.seq\n\n";

$dna_filename = <STDIN>;

chomp $dna_filename;

# Does the file exist?
unless (-e $dna_filename){
	print "File \"$dna_name\" doesnot exist!\n";
	exit;
	}

#Open the file or die
open(DNAFILE, $dna_filename)or die ("I could not open the file $dna_filename!!");

@DNA = <DNAFILE>;

#Open the output file
$psuedo = "psuedo";
open(CHR, ">$psuedo") or die ("I can not create the file!!");

foreach (@DNA){
	if ($_ =~ /^>chr/){
	print CHR "\n$_";
	}
	else{
	chomp $_;
	print CHR "$_";
	}
}
close(CHR);
close(DNAFILE);
exit;
