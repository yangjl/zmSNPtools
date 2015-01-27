=head1 NAME

VertRes::Utils::Cigar - cigar utility functions

=head1 SYNOPSIS

use VertRes::Utils::Cigar;

my $cigar_util = VertRes::Utils::Cigar->new();

# use any of the utility functions described here, eg.


=head1 DESCRIPTION

General utility functions for working on or with cigar output.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Cigar;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Utils::Seq;
use VertRes::Parser::fastq;
use VertRes::Parser::fasta;
use VertRes::Utils::Sam;

use base qw(VertRes::Base);

our $MIN_HIT_LENGTH_454 = 30;

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Cigar->new();
 Function: Create a new VertRes::Utils::Cigar object.
 Returns : VertRes::Utils::Cigar object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 top_10_hits_per_read

 Title   : top_10_hits_per_read
 Usage   : $obj->top_10_hits_per_read($infh, $outfh);
 Function: Filters the cigar output of ssaha2, giving just the top 10 hits per
           read.
 Returns : n/a
 Args    : input (ssaha cigar stream) and output filehandles

=cut

sub top_10_hits_per_read {
    my ($self, $infh, $outfh) = @_;
    
    my $current_read;
    my @current_hits;
    my $finished;
    while (<$infh>) {
        if (/^cigar::/) {
            my @s = split;
            
            $current_read ||= $s[1];
            
            if ($current_read ne $s[1]) {
                for (my $i = 0; $i < 10 && $i < @current_hits; $i++) {
                    print $outfh $current_hits[$i];
                }
                
                $current_read = $s[1];
                undef(@current_hits);
            }
            
            push(@current_hits, $_);
        }
        elsif (/^SSAHA2 finished/) {
            $finished = $_;
        }
    }
    for (my $i = 0; $i < 10 && $i < @current_hits; $i++) {
        print $outfh $current_hits[$i];
    }
    
    if ($finished) {
        print $outfh $finished;
    }
    
    close($infh);
    close($outfh);
    
    return;
}

=head2 cigar_to_sam

 Title   : cigar_to_sam
 Usage   : $obj->cigar_to_sam(\@fastqs, \@cigars, $i_size, $rg, $out_sam);
 Function: Convert ssaha2 cigar output to sam format. Given 2 fastqs and the
           corresponding cigars, pairs up reads correctly in the output sam.
           NB: fastq files must 'match' in that they have the same number of
           sequences in the same order, so mates are the Nth sequence in the
           files.
 Returns : int (the number of sam records written)
 Args    : array ref of fastq files (max 2), array ref of cigar files (max 2,
           the output of running ssaha2 on the corresponding fastqs),
           expected insert size (if supplying 2 fastqs and cigars, else undef),
           read group id (or undef if not known), output file

=cut

sub cigar_to_sam {
    my ($self, $fastqs, $cigars, $insert_size, $read_group, $out_sam) = @_;
    $read_group = $read_group ? 'RG:Z:'.$read_group : '';
    
    my @fastq_parsers;
    my $fasta_input = 0;
    foreach my $fastq (@{$fastqs}) {
        my $parser;
        if ($fastq !~ /q(\.gz)?$/) {
            $fasta_input = 1;
            $parser = VertRes::Parser::fasta->new(file => $fastq);
        }
        else {
            $parser = VertRes::Parser::fastq->new(file => $fastq);
        }
        my $rh = $parser->result_holder;
        push(@fastq_parsers, [$parser, $rh]);
    }
    
    open(my $out_fh, '>', $out_sam) or $self->throw("Cannot create sam file: $!");
    
    my @cigar_hashes;
    foreach my $cigar (@{$cigars}) {
        #*** memory hog...
        push(@cigar_hashes, $self->cigar_by_reads($cigar));
    }
    
    my $seq_util = VertRes::Utils::Seq->new();
    
    my $numReadsWritten = 0;
    my $numReadsPaired = 0;
    my $numReadsDiscarded = 0;
    my $totalReadsHit = 0;
    my ($first_parser, $first_rh) = @{$fastq_parsers[0]};
    while ($first_parser->next_result) {
        my $read_name = $first_rh->[0];
        
        # get the corresponding data from each input fastq
        my (@names, @seqs, @quals);
        foreach my $i (0..$#fastq_parsers) {
            my ($parser, $rh) = @{$fastq_parsers[$i]};
            $parser->next_result unless $i == 0;
            
            my $this_read_name = $read_name;
            if ($rh->[0] ne $read_name) {
                # allow them to differ by the last character, eg. the
                # forward and reverse have been differentiated by name in
                # the fastqs
                if ($read_name =~ /[\/.:]([12ab])$/) {
                    my $this = $1;
                    my $change_to;
                    if ($this eq '1') {
                        $change_to = 2;
                    }
                    elsif ($this eq '2') {
                        $change_to = 1;
                    }
                    elsif ($this eq 'a') {
                        $change_to = 'b';
                    }
                    else {
                        $change_to = 'a';
                    }
                    
                    my $test_read_name = $read_name;
                    $test_read_name =~ s/[12ab]$/$change_to/;
                    
                    if ($rh->[0] eq $test_read_name) {
                        $this_read_name = $test_read_name;
                    }
                }
            }
            
            if ($rh->[0] eq $this_read_name) {
                $names[$i] = $this_read_name;
                $seqs[$i] = $rh->[1];
                $quals[$i] = $rh->[2];
            }
            else {
                $self->throw("read $read_name was in the first fastq, but not in the other at the same position");
            }
        }
        
        my @unmapped;
        
        # get the hits for each fastq
        my @hits;
        foreach my $i (0..$#names) {
            my $hits = $cigar_hashes[$i]{$names[$i]} if $names[$i];
            $hits[$i] = $hits || [];
        }
        
        if ($#{$hits[0]} >= 0 && $#{$hits[1]} >= 0) {
            $numReadsPaired++;
            $totalReadsHit += 2;
            my $mapQual = 0;
            my $flag = 0;
            my $hit1;
            my $hit2;
            
            # from the list of read1 hits and read2 hits - get 1 hit per read        
            my @one_hit_per_read = @{$self->determine_hits($hits[0], $hits[1], $insert_size)};
            my $samCigar1 = $self->ssaha_cigar_to_sam_cigar($one_hit_per_read[0], $seqs[0]);
            my $samCigar2 = $self->ssaha_cigar_to_sam_cigar($one_hit_per_read[1], $seqs[1]);
            my @s1 = split(/\s+/, $one_hit_per_read[0]);
            my @s2 = split(/\s+/, $one_hit_per_read[1]);
            my @flags = @{$self->determine_sam_flag_paired($one_hit_per_read[0], $one_hit_per_read[1], $insert_size)};
            my $mapScore1 = $self->determine_mapping_score($hits[0]);
            my $mapScore2 = $self->determine_mapping_score($hits[1]);
            
            # if paired consistently, double the mapping scores (suggested by RD)
            if ($flags[2] == 1 && $mapScore1 > -1 && $mapScore2 > -1) {
                $mapScore1 = $mapScore1 * 2;
                $mapScore2 = $mapScore2 * 2;
                
                $mapScore1 = 254 if $mapScore1 > 254;
                $mapScore2 = 254 if $mapScore2 > 254;
            }
            
            my $reverse = $s1[4] eq '-';
            my $seq1 = $reverse ? $seq_util->rev_com($seqs[0]) : $seqs[0];
            my $quals1 = $reverse ? reverse($quals[0]) : $quals[0];
            my $seq2 = $reverse ? $seq_util->rev_com($seqs[1]) : $seqs[1];
            my $quals2 = $reverse ? reverse($quals[1]) : $quals[1];
            
            if ($mapScore1 > -1 && $mapScore2 > -1) {
                $self->_print_sam_line($out_fh, \$numReadsWritten, $names[0], $flags[0], $mapScore1, $samCigar1, $seq1, $quals1, $read_group, \@s1, \@s2);
                $self->_print_sam_line($out_fh, \$numReadsWritten, $names[1], $flags[1], $mapScore2, $samCigar2, $seq2, $quals2, $read_group, \@s2, \@s1);
            }
            # just one read aligns
            elsif ($mapScore1 > -1) {
                $flags[0]->{mate_unmapped} = 1;
                delete $flags[0]->{paired_map};
                $self->_print_sam_line($out_fh, \$numReadsWritten, $names[0], $flags[0], $mapScore1, $samCigar1, $seq1, $quals1, $read_group, \@s1);
                $unmapped[1] = $names[1];
            }
            elsif ($mapScore2 > -1) {
                $flags[1]->{mate_unmapped} = 1;
                delete $flags[1]->{paired_map};
                $self->_print_sam_line($out_fh, \$numReadsWritten, $names[1], $flags[1], $mapScore2, $samCigar2, $seq2, $quals2, $read_group, \@s2);
                $unmapped[0] = $names[0];
            }
            # neither read aligns
            else {
                $unmapped[0] = $names[0];
                $unmapped[1] = $names[1];
            }
        }
        elsif ($#{$hits[0]} >= 0 || $#{$hits[1]} >= 0) {
            $totalReadsHit++;
            my ($good_i, $bad_i) = $#{$hits[0]} >= 0 ? (0, 1) : (1, 0);
            
            my $mapScore = $self->determine_mapping_score($hits[$good_i]);
            
            if ($mapScore > -1) {
                my $samCigar = $self->ssaha_cigar_to_sam_cigar($hits[$good_i]->[0], $seqs[$good_i]);
                
                my @s1 = split(/\s+/, $hits[$good_i]->[0]);
                
                my $flag;
                if ($names[0] && $names[1]) {
                    my $flags = $self->determine_sam_flag_paired($hits[0]->[0] || '', $hits[1]->[0] || '');
                    $flag = $flags->[$good_i];
                }
                else {
                    $flag = $self->determine_sam_flag_unpaired($hits[$good_i]->[0]);
                }
                
                my $reverse = $s1[4] eq '-';
                my $seq = $reverse ? $seq_util->rev_com($seqs[$good_i]) : $seqs[$good_i];
                my $quals = $reverse ? reverse($quals[$good_i]) : $quals[$good_i];
                
                $self->_print_sam_line($out_fh, \$numReadsWritten, $names[$good_i], $flag, $mapScore, $samCigar, $seq, $quals, $read_group, \@s1);
            }
            else {
                $unmapped[$good_i] = $names[$good_i];
            }
            
            $unmapped[$bad_i] = $names[$bad_i];
        }
        else {
            $unmapped[0] = $names[0];
            $unmapped[1] = $names[1];
        }
        
        # print out unmapped reads as well
        my @s = (0, 0, 0, 0, 0, '*', 0);
        my $paired_in_tech = $names[0] && $names[1] ? 1 : 0;
        my $flag = {paired_tech => $paired_in_tech,
                    self_unmapped => 1};
        if ($unmapped[0] && $unmapped[1]) {
            $flag->{mate_unmapped} = 1;
        }
        if ($unmapped[0]) {
            $flag->{'1st_in_pair'} = 1 if $paired_in_tech;
            $self->_print_sam_line($out_fh, \$numReadsWritten, $names[0], $flag, 0, '*', $seqs[0], $quals[0], $read_group, \@s);
        }
        if ($unmapped[1]) {
            delete $flag->{'1st_in_pair'};
            $flag->{'2nd_in_pair'} = 1 if $paired_in_tech;
            $self->_print_sam_line($out_fh, \$numReadsWritten, $names[1], $flag, 0, '*', $seqs[1], $quals[1], $read_group, \@s);
        }
    }
    close($out_fh);
    
    return $numReadsWritten;
}

sub _print_sam_line {
    my ($self, $out_fh, $count_ref, $name, $flag_parts, $map_score, $cigar, $seq, $qual, $read_group, $s1, $s2) = @_;
    
    my $sam_util = VertRes::Utils::Sam->new();
    my $flag = $sam_util->calculate_flag(%{$flag_parts});
    
    print $out_fh join("\t", ($name, $flag, $s1->[5], $s1->[6], $map_score, $cigar)), "\t";
    if ($s2) {
        my $equal_field = $s2->[5] eq $s1->[5] ? "=" : $s2->[5];
        my $insert_size = $self->determine_insert_size($s1->[6], $s1->[7], $s2->[6], $s2->[7]);
        print $out_fh join("\t", ($equal_field, $s2->[6], $insert_size)), "\t";
    }
    else {
        print $out_fh join("\t", ('*', 0, 0)), "\t";
    }
    $qual ||= '*';
    print $out_fh join("\t", ($seq, $qual));
    
    if ($read_group) {
        print $out_fh "\t$read_group";
    }
    
    print $out_fh "\n";
    
    $$count_ref++;
}

=head2 cigar_by_reads

 Title   : cigar_by_reads
 Usage   : my $by_reads = $obj->cigar_by_reads($cigar_file);
 Function: Parse a cigar file, group the records according to reads.
 Returns : hash ref where keys are reads and values are array refs
 Args    : cigar file name

=cut

sub cigar_by_reads {
    my ($self, $cigar_file) = @_;
    
    my $cigar_io = VertRes::IO->new(file => $cigar_file);
    my $cigar_fh = $cigar_io->fh;
    
    my %ssaha_read_matches;
    my $ssaha_previous_read = '';
    my $ssaha_current_read_ref = [];
    while (<$cigar_fh>) {
        /^cigar::/ || next;
        chomp;
        
        my (undef, $this_read) = split;
        $ssaha_previous_read ||= $this_read;
        
        if ($ssaha_previous_read ne $this_read) {
            $ssaha_read_matches{$ssaha_previous_read} = [@{$ssaha_current_read_ref}];
            
            $ssaha_current_read_ref = [];
            
            $ssaha_previous_read = $this_read;
        }
        
        push(@{$ssaha_current_read_ref} , $_);
    }
    close($cigar_fh);
    
    # put the entries for the last read in the hash table
    $ssaha_read_matches{$ssaha_previous_read} = $ssaha_current_read_ref;
    
    return \%ssaha_read_matches;
}

=head2 ssaha_cigar_to_sam_cigar

 Title   : ssaha_cigar_to_sam_cigar
 Usage   : my $sam_cigar = $obj->ssaha_cigar_to_sam_cigar($cigar, $seq);
 Function: ??
 Returns : string
 Args    : ??, ??

=cut

sub ssaha_cigar_to_sam_cigar {
    my ($self, $ssaha_cigar, $seq) = @_;
    
    my @s = split(/\s+/, $ssaha_cigar);
    
    # work out the soft clip of the reads (if any)
    my $softClip = '';
    if ($s[4] eq '+' && $s[2] > 1) {
        $softClip = ($s[2] - 1).'S';
    }
    elsif ($s[4] eq '-' && $s[2] < length($seq)) {
        $softClip = (length($seq) - $s[2]).'S';
    }
    
    my $endSoftClip = '';
    if ($s[4] eq '+' && $s[3] < length($seq)) {
        $endSoftClip = (length($seq) - $s[3]).'S';
    }
    elsif ($s[4] eq '-' && $s[3] > 1 ) {
        $endSoftClip = ($s[3] - 1).'S';
    }
    
    # get the cigar string - need to switch the M/I/D for SAM
    my @cigar = splice(@s, 10);
    for (my $i = 0; $i < @cigar; $i += 2) {
        my $t = $cigar[$i];
        $cigar[$i] = $cigar[$i + 1];
        $cigar[$i + 1] = $t;
    }
    
    return $softClip.join( '', @cigar ).$endSoftClip;
}

=head2 determine_hits

 Title   : determine_hits
 Usage   : my $two_hits = $obj->determine_hits($hits1, $hits2, $i_size);
 Function: Determins the best hit per read.
 Returns : array ref of the best 2 hits
 Args    : two array refs of ssaha cigar hits, expected insert size

=cut

sub determine_hits {
    my ($self, $hits1ref, $hits2ref, $expectedInsert) = @_;
    my @hits1 = @{$hits1ref};
    my @hits2 = @{$hits2ref};

    if (@hits1 == 1 && @hits2 == 1) {
        my @hit1Split = split(/\s+/, $hits1[0]);
        my @hit2Split = split(/\s+/, $hits2[0]);
        my $insert = abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3]));
        return [$hits1[0], $hits2[0]];
    }
    elsif (@hits1 == 1 && @hits2 > 1) {
        my @hit1Split = split(/\s+/, $hits1[0]);
        my $hit1Chr = $hit1Split[5];
        
        my $closestOffset = -1;
        my $closestHit = $hits2[0];
        foreach (@hits2) {
            chomp;
            my @hit2Split = split(/\s+/, $_);
            if ($_ =~ /\s+$hit1Chr\s+/ && (abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3])) - $expectedInsert < $closestOffset || $closestOffset == -1)) {
                $closestOffset = abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3])) - $expectedInsert;
                my $closestHit = $_;
            }
        }
        return [$hits1[0], $closestHit];
    }
    elsif (@hits1 > 1 && @hits2 == 1) {
        my @hit2Split = split(/\s+/, $hits2[0]);
        my $hit2Chr = $hit2Split[5];
        
        my $closestOffset = -1;
        my $closestHit = $hits1[0];
        foreach (@hits1) {
            chomp;
            my @hit1Split = split(/\s+/, $_);
            if ($_ =~ /\s+$hit2Chr\s+/ && (abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3] ) ) - $expectedInsert < $closestOffset || $closestOffset == -1)) {
                $closestOffset = abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3])) - $expectedInsert;
                my $closestHit = $_;
            }
        }
        return [$closestHit, $hits2[0]];
    }
    else {
        # pick the best combination of hits on same chr closest to insert size
        my $closestOffset = -1;
        my $hit1 = '';
        my $hit2 = '';
        foreach (@hits1) {
            my @hit1Split = split(/\s+/, $_);
            my $h1 = $_;
            foreach (@hits2) {
                my @hit2Split = split(/\s+/, $_);
                if ($hit1Split[5] eq $hit2Split[5] && (abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3])) - $expectedInsert < $closestOffset || $closestOffset == -1)) {
                    $closestOffset = abs($self->determine_insert_size($hit1Split[2], $hit1Split[3], $hit2Split[2], $hit2Split[3])) - $expectedInsert;
                    $hit1 = $h1;
                    $hit2 = $_;
                }
            }
        }
        
        if (length($hit1) == 0) {
            $hit1 = $hits1[0];
            $hit2 = $hits2[0];
        }
        
        return [$hit1, $hit2];
    }
}

=head2 determine_mapping_score

 Title   : determine_mapping_score
 Usage   : my $score = $obj->determine_mapping_score($hits);
 Function: Determine a mapping score for ssaha cigar lines.
 Returns : int
 Args    : array ref of ssaha hits

=cut

sub determine_mapping_score {
    my ($self, $hits) = @_;
    my @hits = @{$hits};
    
    my @hit1 = split(/\s+/, $hits[0]);
    if (abs($hit1[2] - $hit1[3]) < $MIN_HIT_LENGTH_454) {
        return -1;
    }
    
    if (@hits == 1) {
        return 254;
    }
    else {
        my @hit2 = split(/\s+/, $hits[1]);
        
        # richard durbin suggested this; it is supposed to downweight reads that
        # hit to lots of places due to repetitiveness
        my $d = ($hit1[9] - $hit2[9]) - (2 * scalar(@hits));
        my $total_hits = scalar(@hits);
        
        # we'll behave like MAQ though and give it a mapping quality of at
        # least 0
        if ($d < 0) {
            $d = 0;
        }
        
        return $d < 255 ? $d : 254;
    }
}

=head2 determine_insert_size

 Title   : determine_insert_size
 Usage   : my $size = $obj->determine_insert_size($st1, $end1, $st2, $end2);
 Function: Determine an insert size.
 Returns : int
 Args    : 4 ints (start1, end1, start2, end2)

=cut

sub determine_insert_size {
    my ($self, $st1, $end1, $st2, $end2) = @_;
    
    my $insert = abs($st1 - $st2);
    my $posNegInsert = 0;
    if (abs($st1 - $end2) > $insert) {
        $insert = abs($st1 - $end2);
        $posNegInsert = $st1 - $end2;
    }
    
    if (abs($end1 - $st2) > $insert) {
        $insert = abs($end1 - $st2);
        $posNegInsert = $end1 - $st2;
    }
    
    if (abs($end1 - $end2) > $insert) {
        $insert = abs($end1 - $end2);
        $posNegInsert = $end1 - $end2;
    }
    
    return $posNegInsert;
}

=head2 determine_sam_flag_paired

 Title   : determine_sam_flag_paired
 Usage   : my $data = $obj->determine_sam_flag_paired($cigar1, $cigar2, $i_size);
 Function: See if these cigar lines are paired.
 Returns : array ref of flag hash refs for each input and boolean (true if
           consistently paired for 454)
 Args    : cigar string1, cigar string2, expected insert size

=cut

sub determine_sam_flag_paired {
    my ($self, $cigar1, $cigar2, $expected) = @_;
    
    my $flag1 = {paired_tech => 1, '1st_in_pair' => 1};
    my $flag2 = {paired_tech => 1, '2nd_in_pair' => 1};
    my @s1 = split(/\s+/, $cigar1);
    $s1[4] ||= '';
    my @s2 = split(/\s+/, $cigar2);
    $s2[4] ||= '';
    
    $flag1->{self_reverse} = 1 if $s1[4] eq '-';
    $flag2->{self_reverse} = 1 if $s2[4] eq '-';
    
    my $consistent = 0;
    
    if (length($cigar1) > 0 && length($cigar2) > 0) {
        # if within 3 times the insert size and on same strand (454) - then paired normally
        if ($s1[5] eq $s2[5] && abs(abs($self->determine_insert_size($s1[2], $s1[3], $s2[2], $s2[3])) - $expected) < $expected * 3 && $s1[4] eq $s2[4]) {
            $flag1->{paired_map} = 1;
            $flag2->{paired_map} = 1;
        }
        
        $flag1->{mate_reverse} = 1 if $s2[4] eq '-';
        $flag2->{mate_reverse} = 1 if $s1[4] eq '-';
        
        $consistent = 1;
    }
    elsif (length($cigar1) > 0) {
        $flag1->{mate_unmapped} = 1;
        $flag2->{self_unmapped} = 1;
        $flag2->{self_reverse} = 0;
    }
    elsif (length($cigar2) > 0 ) {
        $flag2->{mate_unmapped} = 1;
        $flag1->{self_unmapped} = 1;
        $flag1->{self_reverse} = 0;
    }
    else {
        $flag1->{self_unmapped} = 1;
        $flag1->{self_reverse} = 0;
        $flag2->{self_unmapped} = 1;
        $flag2->{self_reverse} = 0;
    }
    
    return [$flag1, $flag2, $consistent];
}

=head2 determine_sam_flag_unpaired

 Title   : determine_sam_flag_unpaired
 Usage   : my $flag = $obj->determine_sam_flag_unpaired($cigar);
 Function: Return a sam-style flag indicating if this cigar line is unpaired.
 Returns : hash ref of flags
 Args    : cigar string

=cut

sub determine_sam_flag_unpaired {
    my ($self, $cigar) = @_;
    my @s = split(/\s+/, $cigar);
    
    my $flag = {};
    
    if (length($cigar) > 0) {
        $flag->{self_reverse} = 1 if $s[4] eq '-';
    }
    else {
        $flag->{self_unmapped} = 1;
    }
    
    return $flag;
}

1;
