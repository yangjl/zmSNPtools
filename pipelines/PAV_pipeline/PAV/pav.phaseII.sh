#!/bin/bash

dbpath="/mnt/01/eddyyeh/analyses/databases/blast+"
min_ident="0.90"
min_cvg="0.10"
min_align_length="30"

if [ $# -ne 3 ]
then
	echo ""
	echo "ERROR: Please specify three arguments to execute this script"
	echo "   ARGUMENT 1: Top level directory containing 'phase I' files"
	echo "   ARGUMENT 2: Sample name"
	echo "   ARGUMENT 3: Number of processors"
	echo ""
	exit 0;
else
	dir=$1
	sample=$2
	proc=$3

	# Create symbolig links
	find $dir -name "$sample.*.gte300.fas" -exec ln -s '{}' . \;
	find $dir -name "$sample.*.gte300.*.gff3" -exec ln -s '{}' . \;

	for prefix in `find -name "$sample.*.gte300.fas" | awk '{ sub(/\.gte300\.fas$/, "", $0); print $0; }' | sort -u`
	do
		if [ ! -e "$prefix".01.best_alignments.txt ]
		then
			01-pick_best_alignments.pl -f "$prefix".*.gff3 > "$prefix".01.best_alignments.txt
		fi

		if [ ! -e "$prefix".02.B73_regions_masked.fas ]
		then
			02-mask_B73_regions.pl --breakpoints "$prefix".01.best_alignments.txt --fasta "$prefix".gte300.fas > "$prefix".02.B73_regions_masked.fas
		fi

		if [ ! -e "$prefix".03.iterations.final.fas ]
		then
			03-perform_partial_realignments.pl --fasta "$prefix".02.B73_regions_masked.fas --db "$dbpath" --out "$prefix".03 --ident "$min_ident" --coverage "$min_cvg" --align "$min_align_length" --proc "$proc"
		fi

		if [ ! -e "$prefix".PAV_Contigs.gff3 ]
		then
			04-generate_PAV_gff3.pl --fasta "$prefix".03.iterations.final.fas --source "$prefix" > "$prefix".PAV_Contigs.gff3
		fi

		if [ ! -e "$prefix".PAV_Contigs.B73-Masked.fas ]
		then
			05-extract_sequences.pl --gff3 "$prefix".PAV_Contigs.gff3 --fasta "$prefix".03.iterations.final.fas > "$prefix".PAV_Contigs.B73-Masked.fas
		fi

		if [ ! -e "$prefix".PAV_Contigs.B73-Unmasked.fas ]
		then
			05-extract_sequences.pl --gff3 "$prefix".PAV_Contigs.gff3 --fasta "$prefix".gte300.fas > "$prefix".PAV_Contigs.B73-Unmasked.fas
		fi

		if [ ! -e "$prefix".PAV_Segments.fas ]
		then
			06-extract_PAV_regions.pl --fasta "$prefix".03.iterations.final.fas > "$prefix".PAV_Segments.fas
		fi

		if [ ! -e "$prefix".PAV_Segments.gte300.fas ]
		then
			fastaclean.pl -f "$prefix".PAV_Segments.fas -l 300 > "$prefix".PAV_Segments.gte300.fas
		fi

		# clean up symbolic links
		rm -r "$prefix".gte300.fas
		rm -r "$prefix".gte300.*.gff3
	done
fi
