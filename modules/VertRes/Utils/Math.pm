=head1 NAME

VertRes::Utils::Math - math utility functions

=head1 SYNOPSIS

use VertRes::Utils::Math;

my $math_util = VertRes::Utils::Math->new();

my $median = $math_util->histogram_median({1 => 5, 2 => 10, 3 => 5});

=head1 DESCRIPTION

General utility functions for doing math/stats stuff.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Math;

use strict;
use warnings;

use List::Util qw(sum);
use base qw(VertRes::Base);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Math->new();
 Function: Create a new VertRes::Utils::Math object.
 Returns : VertRes::Utils::Math object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 histogram_median

 Title   : histogram_median
 Usage   : my $median = $obj->histogram_median($histogram_hash_ref);
 Function: Get the median bin of a histogram.
 Returns : int
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut

sub histogram_median {
    my ($self, $hash) = @_;
    
    # find the half-way frequency
    my $total = 0;
    foreach my $freq (values %{$hash}) {
        $total += $freq;
    }
    my $half = sprintf("%0.0f", $total / 2);
    
    # find the corresponding bin
    my $median = 0; 
    my $current = 0;
    foreach my $bin (sort { $a <=> $b } keys %{$hash}) {
        $current += $hash->{$bin};
        if ($current >= $half) {
            $median = $bin;
            last;
        }
    }
    
    return $median;
}

=head2 histogram_quartiles

 Title   : histogram_quartiles
 Usage   : my $median = $obj->histogram_median($histogram_hash_ref);
 Function: Get the quartliles of a histogram.
 Returns : hash with keys q1, q2, q3
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut

sub histogram_quartiles {
    my ($self, $hash) = @_;
    my %quartiles;
    my $total_freq = sum 0, values %{$hash};
    my $freq = 0;
    foreach my $key (sort { $a <=> $b } keys %{$hash}) {
        $freq += $hash->{$key};
        
        if ($freq >= 0.75 * $total_freq ){ 
            $quartiles{q3} = $key;
            $quartiles{q2} = $key unless defined $quartiles{q2};
            $quartiles{q1} = $key unless defined $quartiles{q1};
            last;
        }
        elsif ($freq >= 0.5 * $total_freq and !(defined $quartiles{q2})) {
            $quartiles{q2} = $key;
            $quartiles{q1} = $key unless defined $quartiles{q1};
        }
        elsif ($freq >= 0.25 * $total_freq and !(defined $quartiles{q1})) {
            $quartiles{q1} = $key;
        }
    }
 
    return %quartiles;
}

=head2 histogram_mean

 Title   : histogram_mean
 Usage   : my $mean = $obj->histogram_mean($histogram_hash_ref);
 Function: Get the mean of a histogram
 Returns : float
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut

sub histogram_mean {
    my ($self, $hash) = @_;
    my $total;
    my $count;

    while (my ($k, $v) = each(%$hash)) {
        $total += $k * $v;
        $count += $v;
    }

    return defined $count ? $total / $count : undef;
}

=head2 histogram_stats

 Title   : histogram_stats
 Usage   : my $mean = $obj->histogram_stats($histogram_hash_ref, $mean);
 Function: Get the standard deviation of a histogram
 Returns : hash of statistics, keys are:
           mean, standard_deviation, q1, q2, q3, total (=sum of values)
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins
           Optional mean of histogram, to save calculating it.

=cut

sub histogram_stats {
    my ($self, $hash, $mean) = @_;
    my %stats;
    my $sd = 0;
    my $total = 0;
    $mean = $self->histogram_mean($hash) unless $mean;
    return undef unless defined $mean;
    $stats{mean} = $mean;

    # we use the formula
    # sd^2 = sum( (x_i - mean)^2 ) / n
    while (my ($num, $freq) = each %{$hash}) {
        $sd += $freq * ( ($mean - $num) ** 2);
        $total += $freq;
    }

    $stats{standard_deviation} = ($sd / $total) ** 0.5;
    $stats{total} = $total;
    my %quartiles = $self->histogram_quartiles($hash);
    foreach (qw[q1 q2 q3]){
        $stats{$_} = $quartiles{$_};
    }
    return %stats;
}

=head2 compare_hash_keys

 Title   : compare_hash_keys
 Usage   : my $same = $obj->compare_hash_keys($hash1, $hash2);
 Function: See if two hashes have the same set of keys.
 Returns : boolean
 Args    : 2 hash refs

=cut

sub compare_hash_keys {
    my ($self, $hash1, $hash2) = @_;
    my %hash1 = %{$hash1};
    my %hash2 = %{$hash2};
    
    my $same = 1;
    if (keys %hash1 != keys %hash2) {
        $same = 0;
    }
    else {
        foreach my $key (keys %hash1) {
            last unless exists $hash2{$key};
            delete $hash2{$key};
        }
        if (%hash2) {
            $same = 0;
        }
    }
    
    return $same;
}

1;
