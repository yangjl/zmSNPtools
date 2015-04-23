#!/bin/bash
# Append line to proper spot in file using file split at line to be appended for structure

touch /Users/cjfiscus/Desktop/Testing/pta.txt # create temp file 
cat /Users/cjfiscus/Desktop/Testing/pt1.txt > /Users/cjfiscus/Desktop/Testing/pta.txt # Copy first part of file (ending where line is to be appended) to temp file  
echo "#define MAXPOPS    q      // int number of populations assumed" >> /Users/cjfiscus/Desktop/Testing/pta.txt # Append line to copy of file 
cat /Users/cjfiscus/Desktop/Testing/pta.txt /Users/cjfiscus/Desktop/Testing/pt2.txt > /Users/cjfiscus/Desktop/Testing/mainparams # Merge 2 files together with appended line
rm /Users/cjfiscus/Desktop/Testing/pta.txt # remove temp file created