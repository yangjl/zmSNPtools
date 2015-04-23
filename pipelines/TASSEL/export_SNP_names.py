# Export SNP names from graphing data
# Chris Fiscus
# 3/4/14
# A script that allows you to export a subset of SNPs located in a particular region to a separate file

#####NOTE#####
# In order to use this script, SNP file must contain absolute positions relative to the beginning of the genome
# or relative to the beginning of the chromosome
#
# This script exports SNP name only! To export entire line, see export_selected_regions.py


def main():

    file_name = input("Enter name of input file: ")
    file_name_w = input("Enter name of output file: ")
    print() # blank line for UI
    limit_low = eval(input("Enter lower bound of region to export (bp): "))
    limit_high = eval(input("Enter higher bound of region to export (bp): "))

    
    
    file = open(file_name, "r") # from previous prompt

    x1 = [] # print list

    
    for line in file: # do this for every line in the file

        y = line

        check = y.split() # split line at spaces

        if limit_high >= eval(check[2]) >= limit_low:   # sort lines by chromosome 
            x1.append(str(check[0])) # write line to list 


    file.close() # close input file


    out = open(file_name_w, "w") # open file to be written (change file in " ")
    out.write(' \n'.join(x1)) # write the list to new file
    x1[:] = [] # empty list
    out.close # close out file


 

main()

