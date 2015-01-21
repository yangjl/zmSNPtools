#!/bin/bash

proc=8

if [ $# -eq 0 ]
then
	echo ""
	echo "ERROR: Please specify top level directory containing ABySS assembled contigs"
	echo ""
	exit 0;
else 
	phaseI=`which pav.phaseI.sh`

	if [ `expr length $phaseI` -ne 0 ]
	then
		# Create sumbolic links to fasta *.gte300.fas files
		find $1 -name "*.gte300.fas" -exec ln -s '{}' . \;

		for prefix in `find *.gte300.fas | awk '{ split($0, a, /\./); print a[1]; }' | sort -u`
		do
			echo "sh $phaseI $1 $prefix $proc"
			rm -r "$prefix".*.gte300.fas
		done
	else
		echo ""
		echo "ERROR: Cannot determine the location of script 'pav.phaseI.sh'"
		echo ""
		exit 0;
	fi
fi
