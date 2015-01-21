#!/bin/bash

proc=8

if [ $# -eq 0 ]
then
	echo ""
	echo "ERROR: Please specify the path to the directory containing 'phase II' files"
	echo ""
	exit 0;
else 
	phaseIII=`which pav.phaseIII.sh`

	if [ `expr length $phaseIII` -ne 0 ]
	then
		# Create symbolic links to PAV segments file
		find $1 -name "*.PAV_Segments.gte300.fas" -exec ln -s '{}' . \;

		for prefix in `find *.gte300.fas | awk '{ split($0, a, /\./); print a[1]; }' | sort -u`
		do
			echo "sh $phaseIII $1 $prefix $proc"
			rm -r "$prefix".*.gte300.fas
		done
	else
		echo ""
		echo "ERROR: Cannot determine the location of script 'pav.phaseIII.sh'"
		echo ""
		exit 0;
	fi
fi
