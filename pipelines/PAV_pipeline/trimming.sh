#!/bin/bash

proc=64

for prefix in `find *_1.fastq.gz | awk '{ sub(/_1\.fastq\.gz$/, "", $0); print $0; }' | sort -u`
do
	if [ ! -e "$prefix".trimmed-paired-1.fq.gz ] || [ ! -e "$prefix".trimmed-paired-2.fq.gz ]
	then
		# decompress
		gzip -cd "$prefix"_1.fastq.gz > "$prefix"_1.fastq
		gzip -cd "$prefix"_2.fastq.gz > "$prefix"_2.fastq

		# Trimming
		trim_fastq.pl -i "$prefix"_1.fastq -it fastq-sanger -ot fastq-sanger -o "$prefix".trimmed-1.fq -l "$prefix".trimmed-1.log -p "$proc"
		trim_fastq.pl -i "$prefix"_2.fastq -it fastq-sanger -ot fastq-sanger -o "$prefix".trimmed-2.fq -l "$prefix".trimmed-2.log -p "$proc"
		sort_trimmed.pl -f1 "$prefix".trimmed-1.fq -f2 "$prefix".trimmed-2.fq -p "$prefix".trimmed
		
		rm -r "$prefix"_1.fastq "$prefix"_2.fastq "$prefix".trimmed-1.fq "$prefix".trimmed-2.fq
		screen -d -m gzip -f "$prefix".trimmed-*.log "$prefix".trimmed-*.fq
		sleep 1s
	fi
done
