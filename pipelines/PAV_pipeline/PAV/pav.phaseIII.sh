#!/bin/bash

dbpath="/mnt/01/eddyyeh/analyses/databases/blast+"
ident=0.9
cvg=0.8

if [ $# -ne 3 ]
then
	echo ""
	echo "ERROR: Please specify three arguments to execute this script"
	echo "   ARGUMENT 1: Top level directory containing 'phase II' files"
	echo "   ARGUMENT 2: Sample name"
	echo "   ARGUMENT 3: Number of processors"
	echo ""
	exit 0;
else 
	dir=$1
	sample=$2
	proc=$3

	# Create symbolic links to PAV segments file
	find $dir -name "$sample.*.PAV_Segments.gte300.fas" -exec ln -s '{}' . \;

	for prefix in `find $sample.*.PAV_Segments.gte300.fas | awk '{ sub(/\.fas$/, "", $0); print $0; }' | sort -u`
	do
		if [ ! -e "$prefix".length ]
		then
			sequence_length.pl -f "$prefix".fas > "$prefix".length
		fi
		
		db="fungi-78_Genomes-20130805"
		if [ ! -e "$prefix"."$db".blastn.p"$ident"c"$cvg"it.gff3 ]
		then
			echo "$prefix.fas: Aligning to Fungi Genomes - $db"
			blastn -db "$dbpath"/"$db" -query "$prefix".fas -out "$prefix"."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".length -sl "$dbpath"/"$db".length -blast "$prefix"."$db".blastn.tabular \
				-o "$prefix"."$db".blastn.p"$ident"c"$cvg"it.gff3 -i "$ident" -c "$cvg" -it
		fi

		db="prokaryotes-20130730"
		if [ ! -e "$prefix"."$db".blastn.p"$ident"c"$cvg"it.gff3 ]
		then
			echo "$prefix.fas: Aligning to Prokaryote Genomes - $db"
			blastn -db "$dbpath"/"$db" -query "$prefix".fas -out "$prefix"."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".length -sl "$dbpath"/"$db".length -blast "$prefix"."$db".blastn.tabular \
				-o "$prefix"."$db".blastn.p"$ident"c"$cvg"it.gff3 -i "$ident" -c "$cvg" -it
		fi

		db="homo_sapiens_GRCh37"
		if [ ! -e "$prefix"."$db".blastn.p"$ident"c"$cvg"it.gff3 ]
		then
			echo "$prefix.fas: Aligning to Human genome - $db"
			blastn -db "$dbpath"/"$db" -query "$prefix".fas -out "$prefix"."$db".blastn.tabular \
				-dust yes -num_threads $proc -max_target_seqs 3 -outfmt 6 -task "blastn"
			blastn_tabular2gff3.pl -ql "$prefix".length -sl "$dbpath"/"$db".length -blast "$prefix"."$db".blastn.tabular \
				-o "$prefix"."$db".blastn.p"$ident"c"$cvg"it.gff3 -i "$ident" -c "$cvg" -it
		fi

		echo "$prefix.fas: Removing Contaminants"
		extract_final_PAVs.pl -f "$prefix".fas -g "$prefix".*.p"$ident"c"$cvg"it.gff3 -o "$prefix".final.fas
		sequence_length.pl -f "$prefix".final.fas > "$prefix".final.length

		rm -r "$prefix".fas
	done
fi
