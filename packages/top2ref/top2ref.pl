#!/usr/bin/env perl
# https://github.com/VertebrateResequencing/vr-codebase/blob/master/scripts/topbot-to-fwd-strand
# Author: petr.danecek@sanger
#
# "TOP/BOT" Strand and "A/B" Allele. Illumina SNP Genotyping, technical note.
# Create map file:
#       Input:
#         20  11244   A,C
#         20  11799   A,C
#
#       Output:
#         20  11244  A,C   1   # fwd strand
#         20  11799  A,C  -1   # rev strand, should be flipped
#
# Convert SNP files:
#       Input:
#         20  11244   CA
#         20  11799   AG
#
#       Output:
#         20  11244   CA
#         20  11799   TC
#


use strict;
use warnings;
use Carp;
use FaSlice;

my $opts = parse_params();
if ( exists($$opts{refseq}) ) { create_map($opts); }
elsif ( exists($$opts{map}) ) { convert($opts); }

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    die
        "Usage: topbot-to-fwd-strand [OPTIONS]\n",
        "       cat snps.top | topbot-to-fwd-strand -r ref.fa > snps.map\n",
        "       cat snps.top | topbot-to-fwd-strand -m snps.map > snps.fwd\n",
        "Options:\n",
        "   -m, --map <file>                Output of the previous run with -r.\n",
        "   -r, --refseq <file>             Samtools indexed reference fasta file.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
}


sub parse_params
{
    my $opts = {};
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-m' || $arg eq '--map' ) { $$opts{map}=shift(@ARGV); next }
        if ( $arg eq '-r' || $arg eq '--refseq' ) { $$opts{refseq}=shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{map}) && !exists($$opts{refseq}) ) { error("Missing the -m or -r option.\n") }
    return $opts;
}


sub convert
{
    my ($opts) = @_;
    my %flip = ( A=>'T', C=>'G', G=>'C', T=>'A', N=>'N' );

    open(my $map,'<',$$opts{map}) or error("$$opts{map}: $!");
    while (my $map_line=<$map>)
    {
        my $line = <STDIN>;
        if ( !defined $line ) { error("Too many lines in the map file, starting with $map_line"); }

        chomp($map_line);
        my @map = split(/\t/,$map_line);

        chomp($line);
        my @items = split(/\t/,$line);

        if ( $map[0] ne $items[0] or $map[1] ne $items[1] ) { error("Out of sync: $map[0]:$map[1] vs $items[0]:$items[1]"); }

        my %expected_als = map { $_ => 1 } split(/,/,$map[2]);
        $expected_als{N} = 1;

        my @als = split(//,$items[2]);
        for my $al (@als)
        {
            if ( !exists($expected_als{$al}) ) { error("Alleles not consistent: $map[0]:$map[1] $items[2] vs $map[2]\n"); }
        }

        if ( $map[3] ne '-1' && $map[3] ne '1' ) { error("Neither FWD nor REV? $map_line\n"); }

        my $als = ( $map[3] eq '-1' ) ? $flip{$als[0]}.$flip{$als[1]} : $items[2];
        print "$items[0]\t$items[1]\t$als\n";
    } 
    close($map) or error("close $$opts{map}");
    my $line = <STDIN>;
    if ( defined $line ) { error("Too many lines from STDIN, starting with $line"); }
}


sub create_map
{
    my ($opts) = @_;
    my $refseq = FaSlice->new(file=>$$opts{refseq}, size=>1_000_000);

    while (my $line=<STDIN>)
    {
        # 20  75720   A,C,N
        if ( !($line=~/^(\S+)\t(\d+)\t(\S+)$/) ) { error("Could not parse: $line"); }
        my $chr = $1;
        my $pos = $2;
        my @als = grep { $_ ne 'N' } split(/,/,$3);
        if ( @als != 2 ) { error("Expected two alleles, got [$3] at $chr:$pos"); }
        my $strand = topbot_strand($refseq,$chr,$pos,@als);
        if ( !$strand ) { $strand = 0; }
        
        ### get ref_base
        my $ref_base  = $refseq->get_base($chr,$pos);
        print "$chr\t$pos\t",join(',',@als),"\t$strand", "\t$ref_base\n";
    }
}


sub topbot_strand
{
    my ($refseq,$chr,$pos,@als) = @_;

    my $win = 100;
    my $ref_slice = $refseq->get_slice($chr,$pos-$win,$pos+$win);
    my $ref_base  = $refseq->get_base($chr,$pos);
    if ( substr($ref_slice,$win,1) ne $ref_base ) { error("Bad slice? $chr:$pos $ref_base .. $ref_slice\n"); }

    #  Unambiguous pairs: A/C, A/G, T/C, T/G
    #       - knowledge of reference base at the SNP position is enough to 
    #           determine the strand
    #
    #       TOP      REF   ->  ALLELES   TOP_ON_STRAND
    #       -------------------------------------------
    #       A/C     A or C      A/C         1
    #        "      T or G      T/G        -1
    #       A/G     A or G      A/G         1
    #        "      T or C      T/C        -1
    #
    #
    #  Ambiguous pairs:   A/T, C/G
    #       - sequence walking must be performed (simultaneously upstream and downstream) 
    #           until the first unambiguous pair is encountered. The 5' base determines 
    #           the strand.
    #
    #       TOP    5'REF_BASE   ->  ALLELES   TOP_ON_STRAND
    #       ------------------------------------------------
    #       A/T    A or T             A/T          1
    #        "     C or G             T/A         -1
    #       C/G    A or T             C/G          1
    #        "     C or G             G/C         -1
    #
    my $als = join('',sort @als);

    if ( $als eq 'AC' )
    {
        if ( $ref_base eq 'A' or $ref_base eq 'C' ) { return 1; }
        if ( $ref_base eq 'T' or $ref_base eq 'G' ) { return -1; }
    }
    elsif ( $als eq 'AG' )
    {
        if ( $ref_base eq 'A' or $ref_base eq 'G' ) { return 1; }
        if ( $ref_base eq 'T' or $ref_base eq 'C' ) { return -1; }
    }
    elsif ( $als eq 'AT' or $als eq 'CG' )
    {
        my $i = 1;
        while ($i<$win)
        {
            my $a = substr($ref_slice,$win-$i,1);
            my $b = substr($ref_slice,$win+$i,1);
            my $pair = join('', sort($a,$b));
            if ( $pair eq 'AC' or $pair eq 'AG' or $pair eq 'CT' or $pair eq 'GT' )
            {
                if ( $a eq 'A' or $a eq 'T' ) { return 1; }
                if ( $a eq 'C' or $a eq 'G' ) { return -1; }
            }
            $i++;
        }
    }

    return 0;
}
