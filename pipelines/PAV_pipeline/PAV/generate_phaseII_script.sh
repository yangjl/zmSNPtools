#!/bin/bash

proc=8

if [ $# -eq 0 ]
then
	echo ""
	echo "ERROR: Please specify the path to the directory containing 'phase I' files"
	echo ""
	exit 0;
else 
	phaseII=`which pav.phaseII.sh`

	if [ `expr length $phaseII` -ne 0 ]
	then
		# Create sumbolic links to fasta *.gte300.fas files
		find $1 -name "*.gte300.fas" -exec ln -s '{}' . \;

		for prefix in `find *.gte300.fas | awk '{ split($0, a, /\./); print a[1]; }' | sort -u`
		do
			echo "sh $phaseII $1 $prefix $proc"
			rm -r "$prefix".*.gte300.fas
		done
	else
		echo ""
		echo "ERROR: Cannot determine the location of script 'pav.phaseII.sh'"
		echo ""
		exit 0;
	fi
fi
