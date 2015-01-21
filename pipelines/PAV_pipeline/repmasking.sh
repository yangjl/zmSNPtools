#!/bin/bash

dbpath="/mnt/01/eddyyeh/analyses/databases/gmap-gsnap/"
db="cerRepeats+mtecRepeats"
proc=64

for prefix in `find *.trimmed-*.fq.gz | awk '{ split($0, a, /\./); printf("%s.%s.%s\n", a[1], a[2], a[3]); }' | sort -u`
do
	if [ ! -e "$prefix".trimmed.repmasked-paired-1.fq.gz ]
	then
		# decompress *.trimmed-paired-1.fq.gz
		gzip -cd "$prefix".trimmed-paired-1.fq.gz > "$prefix".trimmed-paired-1.fq
		repeatmaskerFastq.pl -i "$prefix".trimmed-paired-1.fq -o "$prefix".trimmed.repmasked-paired-1.fq \
			-db "$db" -dbpath "$dbpath" -p "$proc"
		rm -r "$prefix".trimmed-paired-1.fq
		screen -d -m gzip "$prefix".trimmed.repmasked-paired-1.fq
	fi

	if [ ! -e "$prefix".trimmed.repmasked-paired-2.fq.gz ]
	then
		# decompress *.trimmed-paired-2.fq.gz
		gzip -cd "$prefix".trimmed-paired-2.fq.gz > "$prefix".trimmed-paired-2.fq
		repeatmaskerFastq.pl -i "$prefix".trimmed-paired-2.fq -o "$prefix".trimmed.repmasked-paired-2.fq \
			-db "$db" -dbpath "$dbpath" -p "$proc"
		rm -r "$prefix".trimmed-paired-2.fq
		screen -d -m gzip "$prefix".trimmed.repmasked-paired-2.fq
	fi

	if [ ! -e "$prefix".trimmed.repmasked-singletons-1.fq.gz ]
	then
		# decompress *.trimmed-singletons-1.fq.gz
		gzip -cd "$prefix".trimmed-singletons-1.fq.gz > "$prefix".trimmed-singletons-1.fq
		repeatmaskerFastq.pl -i "$prefix".trimmed-singletons-1.fq -o "$prefix".trimmed.repmasked-singletons-1.fq \
			-db "$db" -dbpath "$dbpath" -p "$proc"
		rm -r "$prefix".trimmed-singletons-1.fq
		screen -d -m gzip "$prefix".trimmed.repmasked-singletons-1.fq
	fi
	
	if [ ! -e "$prefix".trimmed.repmasked-singletons-2.fq.gz ]
	then
		# decompress *.trimmed-singletons-2.fq.gz
		gzip -cd "$prefix".trimmed-singletons-2.fq.gz > "$prefix".trimmed-singletons-2.fq
		repeatmaskerFastq.pl -i "$prefix".trimmed-singletons-2.fq -o "$prefix".trimmed.repmasked-singletons-2.fq \
			-db "$db" -dbpath "$dbpath" -p "$proc"
		rm -r "$prefix".trimmed-singletons-2.fq
		screen -d -m gzip "$prefix".trimmed.repmasked-singletons-2.fq
	fi

	sleep 1s
done
