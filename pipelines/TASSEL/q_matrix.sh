#!/bin/bash

# This script formats a Q matrix from structure for use with TASSEL
#
# For input file: (must be named input.txt)
# Copy data from structure output in format:

# 1  zheng58    (0)   :  1.000 0.000 0.000 
# 2  zheng22    (0)   :  0.285 0.416 0.299 
# 3     P138    (0)   :  0.115 0.498 0.387 
####################################################################

# Prompt for number of populations  
echo "Enter the number of Populations (K):"
read K
STRING=''

for i in $(seq 1 1 $K);
do		
	STRING="${STRING}	Q$i"
done

cut -c4- input.txt >> output.txt  # removes numbers proceeding the individual name, may need to alter for large number of indiv.

sed "s/(0)//g" output.txt > output2.txt # remove (0) from each line, can be an issue if this is part of individual name

sed "s/://g" output2.txt > output.txt # remove : from each line, can be an issue if this is part of individual name

rm output2.txt # remove temporary file 

echo "<Covariate>" >> output2.txt
echo "<Trait>	$STRING" >> output2.txt 
cat output2.txt output.txt > output3.txt
rm output2.txt 
mv output3.txt output.txt
