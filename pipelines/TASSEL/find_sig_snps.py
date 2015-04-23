# Find Significant SNPs
# Chris Fiscus
# 3/12/14
#
# This script allows you to set a threshold and export SNPs from GWAS results above a given significance threshold
#
# File must be formatted with P value in position 3:
#
# SNP CHR BP P 
# ZMLCHR1P1 1 2802 0.38524991448885704
# ZMLCHR1P2 1 2881 NA
#
# Output will be entire line

def main():

    from math import log 
    
    p = [] # list to be printed

    count = 0

    file = open("mlm_graph_3.txt", "r") # open read file
    
    ##### Locate SNPs within parameters #####
    
    for line in file:

        x = line

        if count == 0: # if first line
            p.append(str(line))
            count = count + 1 # first line done

        else: # all lines after header

            y = x.split() # x split at spaces

            z = y[3] # p value

            if z != "NA": # if p value has a value
            
                z = eval(z) # change to type float

                q = -log(z,10)

                if q >= 6: # significance threshold (10^ -VALUE)
                    p.append(line)


            else: # no p value for this SNP
                pass 

    file.close()

    outfile = open("sig_snps.txt", "w")
    outfile.write(''.join(p))
    outfile.close()

    ##### Write to new file #####


main()
        

        
