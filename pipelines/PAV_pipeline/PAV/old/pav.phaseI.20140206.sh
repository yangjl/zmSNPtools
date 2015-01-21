#!/bin/bash

dbpath="/mnt/01/eddyyeh/analyses/databases/blast+"
ident=0.9
cvg=0.8
proc=64

if [ $# -eq 0 ]
then
	echo ""
	echo "ERROR: Please specify top level directory containing ABySS assembled contigs"
	echo ""
	exit 0;
else 
	# Create sumbolic links to fasta *.gte300.fas files
	find $1 -name "*.gte300.fas" -exec ln -s '{}' . \;

	for prefix in `find *.gte300.fas | awk '{ sub(/\.fas$/, "", $0); print $0; }' | sort -u`
	do
		if [ ! -e "$prefix".length ]
		then
			sequence_length.pl -f "$prefix".fas > "$prefix".length
		fi

		# AGPv3 alignment
		db="AGPv3+mito+chlo"
		if [ ! -e "$prefix".for_AGPv3."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".fas -out "$prefix".for_AGPv3."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".length -sl "$dbpath"/"$db".length -blast "$prefix".for_AGPv3."$db".blastn.tabular \
				-o "$prefix".for_AGPv3."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_AGPv3."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".fas -o "$prefix".for_AGPv2.fas
			sequence_length.pl -f "$prefix".for_AGPv2.fas > "$prefix".for_AGPv2.length
		fi

		# AGPv2 alignment
		db="AGPv2+mito+chlo"
		if [ ! -e "$prefix".for_AGPv2."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".for_AGPv2.fas -out "$prefix".for_AGPv2."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".for_AGPv2.length -sl "$dbpath"/"$db".length -blast "$prefix".for_AGPv2."$db".blastn.tabular \
				-o "$prefix".for_AGPv2."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_AGPv2."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".for_AGPv2.fas -o "$prefix".for_AGPv1.fas
			sequence_length.pl -f "$prefix".for_AGPv1.fas > "$prefix".for_AGPv1.length
		fi

		# AGPv1 alignment
		db="AGPv1+mito+chlo"
		if [ ! -e "$prefix".for_AGPv1."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".for_AGPv1.fas -out "$prefix".for_AGPv1."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".for_AGPv1.length -sl "$dbpath"/"$db".length -blast "$prefix".for_AGPv1."$db".blastn.tabular \
				-o "$prefix".for_AGPv1."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_AGPv1."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".for_AGPv1.fas -o "$prefix".for_WUGSC.fas
			sequence_length.pl -f "$prefix".for_WUGSC.fas > "$prefix".for_WUGSC.length
		fi
		
		# WUGSC BACs alignment
		db="WUGSC_BACs-20130313.pieces"
		if [ ! -e "$prefix".for_WUGSC."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".for_WUGSC.fas -out "$prefix".for_WUGSC."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".for_WUGSC.length -sl "$dbpath"/"$db".length -blast "$prefix".for_WUGSC."$db".blastn.tabular \
				-o "$prefix".for_WUGSC."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_WUGSC."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".for_WUGSC.fas -o "$prefix".for_NatGenetics.fas
			sequence_length.pl -f "$prefix".for_NatGenetics.fas > "$prefix".for_NatGenetics.length
		fi

		# Lai et al. Nature Genetics alignment
		db="Nat_Genet_B73.novel"
		if [ ! -e "$prefix".for_NatGenetics."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".for_NatGenetics.fas -out "$prefix".for_NatGenetics."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".for_NatGenetics.length -sl "$dbpath"/"$db".length -blast "$prefix".for_NatGenetics."$db".blastn.tabular \
				-o "$prefix".for_NatGenetics."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_NatGenetics."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".for_NatGenetics.fas -o "$prefix".for_MAGIv3.fas
			sequence_length.pl -f "$prefix".for_MAGIv3.fas > "$prefix".for_MAGIv3.length
		fi
		
		# MAGIv3.1 alignment
		db="MAGIv3.1.contigs"
		if [ ! -e "$prefix".for_MAGIv3."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".for_MAGIv3.fas -out "$prefix".for_MAGIv3."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".for_MAGIv3.length -sl "$dbpath"/"$db".length -blast "$prefix".for_MAGIv3."$db".blastn.tabular \
				-o "$prefix".for_MAGIv3."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_MAGIv3."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".for_MAGIv3.fas -o "$prefix".for_MAGIv4.fas
			sequence_length.pl -f "$prefix".for_MAGIv4.fas > "$prefix".for_MAGIv4.length
		fi

		# MAGIv4.0 alignment
		db="MAGIv4.0.contigs"
		if [ ! -e "$prefix".for_MAGIv4."$db".blastn.p$ident.c$cvg.it.gff3 ]
		then
			blastn -db "$dbpath"/"$db" -query "$prefix".for_MAGIv4.fas -out "$prefix".for_MAGIv4."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".for_MAGIv4.length -sl "$dbpath"/"$db".length -blast "$prefix".for_MAGIv4."$db".blastn.tabular \
				-o "$prefix".for_MAGIv4."$db".blastn.p$ident.c$cvg.it.gff3 -i $ident -c $cvg -it
			extract_not_aligned.pl -gff3 "$prefix".for_MAGIv4."$db".blastn.p$ident.c$cvg.it.gff3 -fasta "$prefix".for_MAGIv4.fas -o "$prefix".phaseI_not_aligned.fas
			sequence_length.pl -f "$prefix".phaseI_not_aligned.fas > "$prefix".phaseI_not_aligned.length
		fi
	done
fi
