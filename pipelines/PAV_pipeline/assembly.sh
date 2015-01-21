#!/bin/bash

genotype="P131"
output_dir="assembly"
proc=64
aligner="map"
k=85
c=10
s=400
n=10
m=80
p=0.95

# Decompressing files
gzip -cd trimmed.repmasked/P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-1.fq.gz > P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-1.fq
gzip -cd trimmed.repmasked/P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-2.fq.gz > P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-2.fq
gzip -cd trimmed.repmasked/P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq.gz > P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq
gzip -cd trimmed.repmasked/P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq.gz > P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq
gzip -cd trimmed.repmasked/P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-1.fq.gz > P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-1.fq
gzip -cd trimmed.repmasked/P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-2.fq.gz > P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-2.fq
gzip -cd trimmed.repmasked/P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq.gz > P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq
gzip -cd trimmed.repmasked/P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq.gz > P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq

# Paired-End libraries
PE1="P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-1.fq P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-2.fq"
PE2="P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-1.fq P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-2.fq"

# Single-End libraries
SE1="P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq"
SE2="P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq"

# Assembly
time -v -o "$genotype".k"$k"-"$aligner".resources abyss-pe default clean \
	aligner="$aligner" np="$proc" c="$c" s="$s" k="$k" n="$n" j="$proc" m="$m" p="$p" \
	name="$genotype".k"$k"-"$aligner" \
	lib="PE1 PE2" PE1="$PE1" PE2="$PE2" \
	se="$SE1 $SE2"  \
	>> "$genotype".k"$k"-"$aligner".log 2> "$genotype".k"$k"-"$aligner".err

# Rename contigs
rename_abyss_assembly.pl -i "$genotype".k"$k"-"$aligner"-contigs.fa -o "$genotype".k"$k"-"$aligner"-contigs.gte300.fas

# Organization
mkdir -p "$output_dir"/"$genotype"
mv P131.SeqCap-T1_1.k85-map.sh "$output_dir"/"$genotype"
mv "$genotype".k"$k"-"$aligner".* "$output_dir"/"$genotype"
mv "$genotype".k"$k"-"$aligner"-* "$output_dir"/"$genotype"

# Cleanup
rm -r P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-1.fq
rm -r P131.rep1.SeqCap-T1_1.trimmed.repmasked-paired-2.fq
rm -r P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-1.fq
rm -r P131.rep2.SeqCap-T1_1.trimmed.repmasked-paired-2.fq
rm -r P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq
rm -r P131.rep1.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq
rm -r P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-1.fq
rm -r P131.rep2.SeqCap-T1_1.trimmed.repmasked-singletons-2.fq
