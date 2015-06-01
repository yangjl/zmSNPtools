
#!/usr/bin/perl -w# ZSD_segment.pl# Usage: perl fastaread_454reads.pl --i <sequence name list> --d <454 fasta file or trimmed sequence file> #open the ePCR file:print ("what is wrong");open($PCR,$ARGV[0]) or die ("cannot open the ePCR file:");while(<$PCR>){	chomp;	my($name, $chr, $start, $end, $esize, $psize, $diff) =split(/\t\, $_);	print ("$name, $chr, $start \n");}close ($PCR);
