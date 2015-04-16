package SamTools;
use strict;
use Carp;
use File::Spec;
use File::Basename;
use Cwd;
use Digest::MD5  qw(md5 md5_hex md5_base64);
use Utils;

use AssemblyTools;

=pod

=head1 DATA STRUCTURES

$FLAGS

=cut

our $FLAGS = 
{
    'paired_tech'    => 0x0001,
    'read_paired'    => 0x0002,
    'unmapped'       => 0x0004,
    'mate_unmapped'  => 0x0008,
    'reverse_strand' => 0x0010,
    'mate_reverse'   => 0x0020,
    '1st_in_pair'    => 0x0040,
    '2nd_in_pair'    => 0x0080,
    'not_primary'    => 0x0100,
    'failed_qc'      => 0x0200,
    'duplicate'      => 0x0400,
};

sub bam_stat_fofn
{
    croak "Usage: bam_stat bam_fofn\n" unless @_ == 1;
    my $bam_fofn = shift;
    
    croak "Cant find bam_fofn file: $bam_fofn\n" unless -f $bam_fofn;
    
    open( F, $bam_fofn ) or die "cant open bam fofn: $bam_fofn\n";
    while( <F> )
    {
        chomp;
        if( -f $_ )
        {
            my $statF = $_."bamstat";
            my $cmd = qq[bsub -q normal -o $$.bam.o -e $$.bam.e perl -w -e "use SamTools;SamTools::bam_stat( '$_', '$statF');"];
            #print $cmd."\n";
            system( $cmd );
        }
        else
        {
            print "Cant find BAM file: $_\n";
        }
    }
    close( F );
}

sub bam_stat
{
    croak "Usage: bam_stat bam output_file\n" unless @_ == 2;
    my $bam = shift;
    my $output = shift;
    
    my $numReads = 0;
    my $numBases = 0;
    my $insert_mean = 0;
    open( BAM, "samtools view $bam |" ) or die "Cannot open bam file: $bam\n";
    while( <BAM> )
    {
        chomp;
        my @s = split( /\s+/, $_ );
        # only consider reads that actually mapped
        next if ($s[1] & $FLAGS->{unmapped});
        $numBases += length( $s[ 9 ] );
        my $ins = abs( $s[ 8 ] );
        $insert_mean = ( ( $insert_mean * $numReads ) + $ins ) / ( $numReads + 1 ) unless $ins > 10000 && $ins != 0;
        $numReads ++;
    }
    close( BAM );
    
    #do a second pass to get the std deviations (uses too much memory if store all inserts in memory)
    my $sumSqVariances = 0;
    open( BAM, "samtools view $bam |" ) or die "Cannot open bam file: $bam\n";
    while( <BAM> )
    {
        chomp;
        my @s = split( /\s+/, $_ );
        my $ins = abs( $s[ 8 ] );
        $sumSqVariances += ( ( $insert_mean - $ins ) * ( $insert_mean - $ins ) ) unless $ins > 10000 && $ins != 0;
    }
    close( BAM );
    
    my $std = $sumSqVariances ? sqrt( ( $sumSqVariances / $numReads ) ) : 0;
    
    open( O, ">$output" ) or die "Cannot create output: $!\n";
    print O "num_bases:$numBases\n";
    print O "num_reads:$numReads\n";
    print O "avg_insert:$insert_mean\n";
    print O "std_dev:$std\n";
    close( O );
}

sub parse_bam_header_line
{
    my ($line) = @_;
    my $out = {};

    # Header line:
    #   @RG     ID:ERR001773    PL:ILLUMINA     PU:IL23_337_6   LB:g1k-sc-NA12878-CEU-2 PI:200  SM:NA12878      CN:SC
    my @items = split /\t/, $line;
    shift(@items);
    for my $pair (@items)
    {
        my ($key,$value) = split /:/,$pair;
        $$out{$key} = $value;
    }
    return $out;
}

sub parse_bam_line
{
    my ($line) = @_;
    my $out = {};

    chomp($line);

    #   IL14_1902:3:83:1158:1446        89      1       23069154        37      54M     =       23069154        0       TGCAC ... RG:Z:ERR001720
    my @items = split /\t/, $line;

    $$out{'flag'}   = $items[1];
    $$out{'chrom'}  = $items[2];
    $$out{'pos'}    = $items[3];
    $$out{'cigar'}  = $items[5];
    $$out{'isize'}  = $items[8];
    $$out{'seq'}    = $items[9];
    $$out{'qual'}   = $items[10];

    my $nitems = @items;
    for (my $i=11; $i<$nitems; $i++)
    {
        my ($key,@vals) = split /:/,$items[$i];

        # e.g. RG:Z:ERR001720
        if ( $key eq 'RG' )
        {
            $$out{'RG'} = $vals[1];
        }

        # e.g. NM:i:0
        elsif( $key eq 'NM' )
        {
            $$out{'NM'} = $vals[1];
        }

        # e.g. MF:i:64
        elsif( $key eq 'MF' )
        {
            $$out{'MF'} = $vals[1];
        }
    }
    return $out;
}

sub soft_clip
{
    my ($threshold,$quals,$len) = @_;

    my $sum = 0;
    my $max = 0;
    my $max_l = $len-1;
    my $l = $max_l;
    for (my $l=$max_l; $l>=35-1; $l--)
    {
        $sum += $threshold - (ord(substr($quals,$l,1)) - 33);
        if ( $sum<0 ) { last; }
        if ( $sum>$max )
        {
            $max = $sum;
            $max_l = $l;
        }
    }
    return ($max_l + 1);
}


=head2 collect_detailed_bam_stats

        Description: Reads the given bam file and runs detailed statistics, such as 
                        histograms of insert size; QC content of mapped and both mapped
                        and unmapped sequences; reads distribution with respect to
                        chromosomes; duplication rate.
        Arg [1]    : The .bam file (sorted and indexed) or a filehandle.
        Arg [2]    : The .fai file (to determine chromosome lengths). (Can be NULL with do_chrm=>0.)
        Arg [3]    : Options (hash) [Optional]
        Options    : (see the code for default values)
                        do_chrm          .. should collect the chromosome distrib. stats.
                        do_clipped       .. collect the number of bases left after BWA clipping, use the given threshold
                        do_gc_content    .. should we collect the gc_content (default is 1)
                        do_rmdup         .. default is 1 (calculate the rmdup)

                        insert_size_bin  .. the length of the distribution intervals for the insert size frequencies
                        gc_content_bin   
                        main_bulk        .. the minimum x-value which covers at least 95% hits in the histograms

        Returntype : Hash with the following entries, each statistics type is a hash as well. See also Graphs::plot_stats.
                     The stats are collected for all reads ('total') and for each reading group (RG) separately.
                        insert_size =>
                            yvals      .. array of insert size frequencies 
                            xvals      .. 
                            max => x   .. maximum values
                            max => y
                            average    .. average value (estimated from histogram, the binsize may influence the accuracy)
                            std_dev    .. standard deviation
                        gc_content_forward  .. gc content of sequences with the 0x0040 flag set
                        gc_content_reverse  .. gc content of seqs with the 0x0080 flag set
                        reads_chrm_distrib  .. read distribution with respect to chromosomes
                        reads_total
                        reads_paired
                        reads_mapped        .. paired + unpaired
                        reads_unpaired
                        reads_unmapped
                        bases_total
                        bases_mapped_read   .. number of mapped reads * read length
                        bases_mapped_cigar  .. number of M+I in cigars
                        duplication
                        num_mismatches (if NM fields defined in bam file, otherwise 0)
=cut

sub collect_detailed_bam_stats
{
    my ($bam_file,$fai_file,$options) = @_;
    if ( !$bam_file ) { Utils::error("Expected .bam file as a parameter.\n") }
    
    $options = {} unless $options;
    my $insert_size_bin = exists($$options{'insert_size_bin'}) ? $$options{'insert_size_bin'} : 1;
    my $gc_content_bin  = exists($$options{'gc_content_bin'}) ? $$options{'gc_content_bin'} : 1;

    my $do_clipped   = exists($$options{'do_clipped'}) ?  $$options{'do_clipped'} : 0;
    my $do_chrm      = exists($$options{'do_chrm'}) ?  $$options{'do_chrm'} : '^(?:\d+|X|Y)$';
    my $chrom_regex  = $do_chrm ? qr/$do_chrm/ : undef;
    my $do_gc        = exists($$options{'do_gc_content'}) ? $$options{'do_gc_content'} : 1;
    my $do_rmdup     = exists($$options{'do_rmdup'}) ? $$options{'do_rmdup'} : 1;
    my $chrm_lengths = $do_chrm ? Utils::fai_chromosome_lengths($fai_file) : {};
    my $main_bulk    = exists($$options{main_bulk}) ? $$options{main_bulk}/100. : 0.95;


    # Use hashes, not arrays - the data can be broken and we might end up allocating insanely big arrays 

    #   reads_unmapped      ..   That is, not aligned to the ref sequence (unmapped flag)
    #   reads_paired        ..      both mates are mapped and form a pair (read_paired flag)
    #   reads_unpaired      ..      both mates mapped, but do not form a pair (neither unmapped nor read_paired flag set)
    #   bases_total         ..   The total number of bases determined as \sum_{seq} length(seq)
    #   bases_mapped_cigar  ..      number of bases with 'M+I' in cigar
    #   bases_mapped_read   ..      number of mapped reads * read length

    my $out_stats = {};

    # Collect the statistics - always collect the total statistics for all lanes ('total') and
    #   when the RG information is present in the header, collect also statistics individually
    #   for each ID (@RG ID:xyz). The statistics are collected into out_stats hash, individual
    #   keys are 'total' and IDs.
    #
    my $i=0;
    my $fh;
    if ( ref($bam_file) eq 'GLOB' ) { $fh=$bam_file; }
    else { open($fh, "samtools view -h $bam_file |") or Utils::error("samtools view -h $bam_file |: $!"); }
    while (my $line=<$fh>)
    {
        # Header line:
        #   @RG     ID:ERR001773    PL:ILLUMINA     PU:IL23_337_6   LB:g1k-sc-NA12878-CEU-2 PI:200  SM:NA12878      CN:SC
        #
        # Data line:
        #   IL14_1902:3:83:1158:1446        89      1       23069154        37      54M     =       23069154        0       TGCAC ... RG:Z:ERR001720
        #
        if ( $line=~/^\@/ )
        {
            if ( !($line=~/^\@RG/) ) { next }

            my $header = parse_bam_header_line($line);
            if ( !exists($$header{'ID'}) ) { Utils::error("No ID in the header line? $line\n"); }

            $$out_stats{$$header{'ID'}}{'header'} = $header;
            next;
        }
		
        # The @stats array is a convenient way how to reuse the same code for the total and individual
        #   statistics - we will always add the same numbers to 'total' and the ID.
        #
        my @stats = ('total');
        my $data  = parse_bam_line($line);
        if ( exists($$data{'RG'}) ) { push @stats, $$data{'RG'}; }

        my $flag   = $$data{'flag'};
        my $chrom  = $$data{'chrom'};
        my $pos    = $$data{'pos'};
        my $cigar  = $$data{'cigar'};
        my $isize  = $$data{'isize'};
        my $seq    = $$data{'seq'};
        
        my $seq_len  = length($seq);
        for my $stat (@stats)
        {
            $$out_stats{$stat}{'reads_total'}++;
            $$out_stats{$stat}{'bases_total'} += $seq_len;
        }

        my $soft_clip = 0;
        my $nmapped = 0;
        my $cigar_info;
        if ( $cigar ne '*' )
        {
            $cigar_info = cigar_stats($cigar);
            if ( exists($$cigar_info{'M'}) ) { $nmapped += $$cigar_info{'M'}; }
            if ( exists($$cigar_info{'I'}) ) { $nmapped += $$cigar_info{'I'}; }
            $soft_clip += $nmapped;
        }
        elsif ( $do_clipped )
        {
            $soft_clip += soft_clip($do_clipped,$$data{qual},$seq_len);
        }
        for my $stat (@stats) { $$out_stats{$stat}{clip_bases} += $soft_clip; } 

        if ( !($flag & $$FLAGS{'unmapped'}) )
        {
            # Stats which make sense only for mapped reads.
            my $mismatch = exists($$data{'NM'}) ? $$data{'NM'} : 0;
            for my $stat (@stats) 
            { 
                # There are two ways how to determine the number of mapped reads:
                #   reads_maped * read_length OR count M+I in cigars. 
                # Let's do both.
                $$out_stats{$stat}{'bases_mapped_read'}  += $seq_len; 
                $$out_stats{$stat}{'bases_mapped_cigar'} += $nmapped; 
                $$out_stats{$stat}{'num_mismatches'} += $mismatch;
            }

            # Chromosome distribution
            #
            if ( $do_chrm && $chrom=~$chrom_regex ) 
            {
                for my $stat (@stats) { $$out_stats{$stat}{'chrm_distrib_freqs'}{$chrom}++; }
            }
        }

        # For long inserts, the libraries are prepared differently and the pairs
        #   point in other direction. bwa does not set the paired flag for them.
        #   Similar thing is true also for 454 reads. Therefore, do the insert
        #   size stats for all unmapped reads above.
        #
        # This is the original version which does not work with long insert sizes:
        #   my $paired = ($flag & $$FLAGS{'read_paired'}) && ($flag & $$FLAGS{'paired_tech'});
        #   if ( $paired )
        #
        # This version gave about half of the paired reads only, the isize can be negative:
        #   my $paired = ($flag & $$FLAGS{'paired_tech'}) && !($flag & $$FLAGS{'unmapped'});
        #   if ( $paired && $isize>0 )

        my $paired = ($flag & $$FLAGS{'paired_tech'}) && !($flag & $$FLAGS{'unmapped'}) && !($flag & $$FLAGS{'mate_unmapped'});
        if ( $paired && $isize!=0 )
        {
            for my $stat (@stats) { $$out_stats{$stat}{'reads_paired'}++; }

            # Insert Size Frequencies
            my $bin = abs(int($isize/$insert_size_bin));
            for my $stat (@stats) { $$out_stats{$stat}{'insert_size_freqs'}{$bin}++; }
        }
        elsif ( $flag & $$FLAGS{'unmapped'} ) 
        { 
            for my $stat (@stats) { $$out_stats{$stat}{'reads_unmapped'}++;  }
        }
        else 
        { 
            for my $stat (@stats) { $$out_stats{$stat}{'reads_unpaired'}++; }
        }

        # GC Content Frequencies - collect stats for both pairs separately. Dont' attempt
        #   to calculate first % GC content - this would occasionally produce spikes. We
        #   work with discrete data. Do the x-axis scaling at the end.
        if ( $do_gc )
        {
            my $gc_count = 0;
            for (my $ipos=0; $ipos<$seq_len; $ipos++)
            {
                my $nuc  = substr($seq, $ipos, 1);
                if ( $nuc eq 'g' || $nuc eq 'G' || $nuc eq 'c' || $nuc eq 'C' ) { $gc_count++; }
            }
            if ( $gc_count )
            {
                if ( $flag & $$FLAGS{'1st_in_pair'} )
                {
                    for my $stat (@stats) 
                    { 
                        $$out_stats{$stat}{'gc_content_fwd_freqs'}{$gc_count}++; 
                        $$out_stats{$stat}{'gc_content_fwd_reads'}++; 
                        $$out_stats{$stat}{'gc_content_fwd_bases'}+= $seq_len; 
                    }
                }
                elsif ( $flag & $$FLAGS{'2nd_in_pair'} )
                {
                    for my $stat (@stats) 
                    { 
                        $$out_stats{$stat}{'gc_content_rev_freqs'}{$gc_count}++; 
                        $$out_stats{$stat}{'gc_content_rev_reads'}++; 
                        $$out_stats{$stat}{'gc_content_rev_bases'}+= $seq_len; 
                    }
                }
                else
                { 
                    # Either it is a non-paired-read technology, or the 1st_in_pair and
                    #   and 2nd_in_pair flags got lost in the process. In that case, 
                    #   both 1st_in_pair and 2nd_in_pair should be set to zero.
                    for my $stat (@stats) 
                    { 
                        $$out_stats{$stat}{'gc_content_single_freqs'}{$gc_count}++; 
                        $$out_stats{$stat}{'gc_content_single_reads'}++; 
                        $$out_stats{$stat}{'gc_content_single_bases'}+= $seq_len; 
                    }
                }
            }
        }

        #if ( $i++>5000 ) { last }
    }
    close $fh;

    # If some type of reads were not present, set the values to 0.
    for my $stat (keys %$out_stats) 
    { 
        $$out_stats{$stat}{'reads_unmapped'} = 0 unless exists($$out_stats{$stat}{'reads_unmapped'});
        $$out_stats{$stat}{'reads_unpaired'} = 0 unless exists($$out_stats{$stat}{'reads_unpaired'});
        $$out_stats{$stat}{'reads_paired'}   = 0 unless exists($$out_stats{$stat}{'reads_paired'});
        $$out_stats{$stat}{'reads_total'}    = 0 unless exists($$out_stats{$stat}{'reads_total'});
        $$out_stats{$stat}{'bases_total'}    = 0 unless exists($$out_stats{$stat}{'bases_total'});
        $$out_stats{$stat}{'bases_mapped_read'}  = 0 unless exists($$out_stats{$stat}{'bases_mapped_read'});
        $$out_stats{$stat}{'bases_mapped_cigar'} = 0 unless exists($$out_stats{$stat}{'bases_mapped_cigar'});

        $$out_stats{$stat}{'reads_mapped'} = $$out_stats{$stat}{'reads_total'} - $$out_stats{$stat}{'reads_unmapped'};
    }

    # Find out the duplication rate. For this, run rmdup.
    if ( $do_rmdup )
    {
        # Run the same routine recursively to get the same stats from the rmdupped file.
        my $rmdup_opts = { do_rmdup=>0, do_chrm=>0, do_gc=>0, do_clipped=>0 };
        open(my $rmdup_fh,"samtools rmdup $bam_file - | samtools view - |") or Utils::error("samtools rmdup $bam_file - | samtools view - |: $!");
        my $rmdup_stats = collect_detailed_bam_stats($rmdup_fh,0,$rmdup_opts);
        close($rmdup_fh);

        for my $stat (keys %$rmdup_stats)
        {
            $$out_stats{$stat}{rmdup_bases_mapped_cigar} = $$rmdup_stats{$stat}{bases_mapped_cigar};
            $$out_stats{$stat}{rmdup_bases_total}        = $$rmdup_stats{$stat}{bases_total};
            $$out_stats{$stat}{rmdup_reads_total}        = $$rmdup_stats{$stat}{reads_total};
            $$out_stats{$stat}{rmdup_reads_mapped}       = $$rmdup_stats{$stat}{reads_mapped};
        }
    }

    # Now process the results. The out_stats hash now contains the total statistics (the key 'total')
    #   and possibly also separate statistics for individual libraries (other keys, e.g. 'ERR001773').
    #   Calculate some numbers (error_rate and duplication) and include bin_size for the histograms
    #   (insert_size and gc_content).
    #
    for my $stat (keys %$out_stats)
    {
        $$out_stats{$stat}{error_rate} = $$out_stats{$stat}{'num_mismatches'} ? 
                                                $$out_stats{$stat}{'num_mismatches'}/$$out_stats{$stat}{'bases_mapped_cigar'} : 0;

        if ( exists($$out_stats{$stat}{'rmdup_bases_mapped_cigar'}) )
        {
            $$out_stats{$stat}{'duplication'} = 1 - $$out_stats{$stat}{'rmdup_reads_mapped'} / $$out_stats{$stat}{'reads_mapped'};
        }

        if ( exists($$out_stats{$stat}{insert_size_freqs}) )
        {
            $$out_stats{$stat}{insert_size} =
            {
                'data'       => $$out_stats{$stat}{'insert_size_freqs'},
                'bin_size'   => $insert_size_bin,
            };
            delete($$out_stats{$stat}{insert_size_freqs});
        }

        if ( exists($$out_stats{$stat}{'gc_content_fwd_freqs'}) )
        {
            $$out_stats{$stat}{'gc_content_forward'} = 
            {
                data        => $$out_stats{$stat}{'gc_content_fwd_freqs'},
                scale_x     => 100.*$$out_stats{$stat}{gc_content_fwd_reads}/$$out_stats{$stat}{gc_content_fwd_bases},
                bin_size    => 1,
            };
            delete($$out_stats{$stat}{'gc_content_fwd_freqs'});
        }
        if ( exists($$out_stats{$stat}{'gc_content_rev_freqs'}) )
        {
            $$out_stats{$stat}{'gc_content_reverse'} =
            {
                data        => $$out_stats{$stat}{'gc_content_rev_freqs'},
                scale_x     => 100.*$$out_stats{$stat}{gc_content_rev_reads}/$$out_stats{$stat}{gc_content_rev_bases},
                bin_size    => 1,
            };
            delete($$out_stats{$stat}{'gc_content_rev_freqs'});
        }
        if ( exists($$out_stats{$stat}{'gc_content_single_freqs'}) )
        {
            $$out_stats{$stat}{'gc_content_single'} =
            {
                data        => $$out_stats{$stat}{'gc_content_single_freqs'},
                scale_x     => 100.*$$out_stats{$stat}{gc_content_single_reads}/$$out_stats{$stat}{gc_content_single_bases},
                bin_size    => 1,
            };
            delete($$out_stats{$stat}{'gc_content_single_freqs'});
        }
    }

    # Convert the hashes into arrays and find extreme values - this is for convenience only,
    #   TrackQC uses this. Although the code is lengthy, the histograms are small and take
    #   no time to process.
    #
    for my $stat_name (keys %$out_stats)
    {
        my $stat = $$out_stats{$stat_name};

        for my $key (keys %$stat)
        {
            if ( ref($$stat{$key}) ne 'HASH' ) { next }
            if ( !exists($$stat{$key}{'data'}) ) { next }
            if ( !exists($$stat{$key}{'bin_size'}) ) { next }

            my @yvals = ();
            my @xvals = ();
            my ($ymax,$xmax);

            my $avg  = 0;
            my $navg = 0;

            my $data = $$stat{$key}{'data'};
            for my $ibin (sort {$a<=>$b} keys %$data)
            {
                my $bin = $ibin * $$stat{$key}{'bin_size'};
                if ( exists($$stat{$key}{scale_x}) ) { $bin *= $$stat{$key}{scale_x}; }
                push @xvals, $bin;

                my $yval = $$data{$ibin};
                push @yvals, $yval;

                $avg  += $yval*$bin;    # yval is the count and bin the value
                $navg += $yval;

                if ( !defined $ymax || $yval>$ymax ) { $xmax=$bin; $ymax=$yval }
            }
            $avg = $navg ? $avg/$navg : 0;

            if ( !@xvals ) { @xvals=(0); }  # Yes, this can happen,e.g. AKR_J_SLX_200_NOPCR_1/1902_3
            if ( !@yvals ) { @yvals=(0); }
            
            my $dev  = 0;
            my $bulk = 0;
            for (my $i=0; $i<scalar @xvals; $i++)
            {
                $dev  += $yvals[$i]*($xvals[$i]-$avg)**2;    # yval is the count and xval the value
                $bulk += $yvals[$i];
                if ( $bulk/$navg>$main_bulk ) { $bulk=$xvals[$i]; last; }
            }
            $dev = $navg ? $dev/$navg : 0;

            $$stat{$key}{'xvals'} = \@xvals;
            $$stat{$key}{'yvals'} = \@yvals;
            $$stat{$key}{'max'}{'x'} = defined $xmax ? $xmax : 0;
            $$stat{$key}{'max'}{'y'} = defined $ymax ? $ymax : 0;
            $$stat{$key}{'average'}  = $avg;
            $$stat{$key}{'std_dev'}  = sqrt($dev);
            $$stat{$key}{main_bulk}  = $bulk;
        }

        # Chromosome distribution histograms (reads_chrm_distrib) are different - xvalues are not numeric
        #
        if ( exists($$stat{'chrm_distrib_freqs'}) )
        {
            my @yvals = ();
            my @xvals = ();
            for my $key (sort Utils::cmp_mixed keys %{$$stat{'chrm_distrib_freqs'}})
            {
                if ( !exists($$chrm_lengths{$key}) ) { Utils::error("The chromosome \"$key\" not in $fai_file.\n") }
                push @yvals, $$stat{'chrm_distrib_freqs'}{$key} / $$chrm_lengths{$key};
                push @xvals, $key;
            }   
            $$stat{'reads_chrm_distrib'}{'yvals'} = scalar @yvals ? \@yvals : [0];
            $$stat{'reads_chrm_distrib'}{'xvals'} = scalar @xvals ? \@xvals : [0];
            $$stat{'reads_chrm_distrib'}{'scaled_dev'} = scaled_chrm_dev(\@xvals,\@yvals);
        }
    }
    
    return $out_stats;
}


# Calculates some sort of standard devitation, except that the
#   the difference from mean is scaled by mean (so that the formula
#   works for samples of different size) and the instead of mean
#   the maximum is used (we want to detect one oversampled chromosome): 
#       sqrt[(1/N)*sum_i((y_i-max)/max)**2]
#
sub scaled_chrm_dev
{
    my ($xvals,$yvals) = @_;
    my $max   = 0;
    my $ndata = 0;
    for (my $i=0; $i<scalar @$xvals; $i++)
    {
        # We do not know if both chromosome X and Y should be counted - we don't know the sex
        if ( !($$xvals[$i]=~/^\d+$/) ) { next } 

        if ( !defined($max) || $max<$$yvals[$i] ) { $max=$$yvals[$i]; }
        $ndata++;
    }
    if ( !$ndata ) { return 0 }

    my $stddev = 0;
    for (my $i=0; $i<scalar @$xvals; $i++)
    {
        if ( !($$xvals[$i]=~/^\d+$/) ) { next } 

        $stddev += (($$yvals[$i] - $max)/$max)**2;
    }
    return sqrt($stddev/$ndata);
}


=head2 cigar_stats

        Description: Parse the SAM cigar and return some info.
        Arg [1]    : The cigar (e.g. 36M, 5M1I30M, 21M15S, etc.)
        Returntype : The hash e.g. { 'M'=>35, 'I'=>1 } 

=cut

sub cigar_stats
{
    my ($cigar) = @_;

    my $stats = {};
    while ($cigar)
    {
        if ( !($cigar=~/^(\d+)(\D)/) ) { last }
        $$stats{$2} += $1;
        $cigar = $';
    }

    return $stats;
}


=head2 print_flags

        Description: For debugging purposes, prints flags for all sequences in human readable form.
        Arg [1]    : The .bam file (sorted and indexed).
        Returntype : None

=cut

sub print_flags
{
    my ($bam_file,$options) = @_;

    if ( !$options ) { $options='' }

    my $fh = \*STDIN;
    if ( $bam_file ) { open($fh, "samtools view $bam_file $options |") or Utils::error("samtools view $bam_file $options |: $!"); }
    while (my $line=<$fh>)
    {
        # IL14_1902:3:83:1158:1446        89      1       23069154        37      54M     =       23069154        0       TGCAC
        my @items = split /\t/, $line;
        my $flag  = $items[1];
        my $chrom = $items[2];
        my $pos   = $items[3];
        my $qual  = $items[4];
        my $isize = $items[8];
        my $seq   = $items[9];

        print $line;
        print "\t insert_size=$isize flag=$flag qual=$qual\n\t --\n";
        print debug_flag($flag);
        print "\n";
    }
    close $fh;
    return;
}

sub debug_flag
{
    my ($flag) = @_;

    my $out = '';
    for my $key (sort keys %$FLAGS)
    {
        if ( $flag & $$FLAGS{$key} ) { $out .= "\t $key\n" }
    }
    return $out;
}


=head2 determineUnmappedFlag

    Arg [1]    : 1/0 flag to say whether the lane is a paired lane
    Arg [2]    : -1/0/1 flag of whether the read is the first read of a pair
    Example    : determineUnmappedFlag( 0, 1 );
    Description: works out the SAM flag value for an unmapped read
    Returntype : int
=cut
sub determineUnmappedFlag
{
    croak "Usage: determineUnmappedFlag paired(0/1) read1(-1/0/1)\n" unless @_ == 2;
    
    my $paired = shift;
    my $read1 = shift;
    
    croak "Paired flag must be 0 or 1\n" unless $paired == 0 || $paired == 1;
    
    croak "Read1 flag must be -1, 0 or 1\n" unless $read1 == 0 || $read1 == 1 || $read1 == -1;
    
    my $total = $paired == 1 ? hex( '0x0001' ) : 0;
    $total += hex( '0x0004' );
    $total += hex( '0x0008' );
    
    if( $paired == 1 )
    {
        #if the read is a fragment from an paired run - then it wont be either end of a pair (e.g. 454)
        if( $read1 != -1 )
        {
            $total += $read1 == 1 ? hex( '0x0040' ) : hex( '0x0080' );
        }
    }
    
    return $total;
}

=head2 pileup2Intervals

    Arg [1]    : samtools pileup file
    Arg [2]    : sam header file
    Arg [3]    : output intervals file
    Example    : pileup2Intervals
    Description: creates an intervals file for the broad recalibrator which excludes all positions +/- 20bp around all indels
    Returntype : none
=cut

sub pileup2Intervals
{
    croak "Usage: pileup2Intervals samtools_pileup_file sam_header_file output_file\n" unless @_ == 3;
    
    my $pileup = shift;
    my $sam_header = shift;
    my $output = shift;
    
    croak "Cant find pileup file: $pileup" unless -f $pileup;
    croak "Cant find sam header file: $sam_header" unless -f $sam_header;
    
    my $pfh;
    if( $pileup =~ /\.gz$/ )
    {
        open( $pfh, "gunzip -c $pileup |" ) or die "Cannot open pileup file: $!\n";
    }
    else
    {
        open( $pfh, $pileup ) or die "Cannot open pileup file: $!\n";
    }
    
    system( qq[cat $sam_header > $output] );
    
    open( my $out, ">>$output" ) or die "Cannot create output file: $!\n";
    my $currentStartPos = 1;
    my $c = 0;
    
    my $currentChr = 'none';
    while( <$pfh> )
    {
        chomp;
        if( $_ =~ /^.+\t\d+\t\*\t.*/ )
        {
            my @s = split( /\t/, $_ );
            
            if( $currentChr eq 'none' || $currentChr ne $s[ 0 ] )
            {
                $currentChr = $s[ 0 ];
                $currentStartPos = 1;
            }
            
            my $stop = $s[ 1 ] - 20;
            
            if( $stop > $currentStartPos )
            {
                print $out qq/$s[ 0 ]\t$currentStartPos\t$stop\t+\ttarget_$c\n/;
                $currentStartPos = $s[ 1 ] + 20;
                $c ++;
            }
        }
    }
    close( $pfh );
    close( $out );
}

1;
