# Export SNP names for multiple genomic regions
# Chris Fiscus
# 3/7/14
#
# NOTE: Process .csv file with $sed "s/^M/\n/g" ./filename.csv > output.txt before using this script. Use output.txt as first input. 
#
# This script takes a .txt file formatted:
#
#Chr,Start,Stop,Midpoint,Cluster Size,Distance to Next
#Chr01,134858461,134868286,134863373.5,9825,29686
#Chr07,22949640,22959358,22954499,9718,14062
#
#and reads which regions to use from a particular chromosome (designated in line 36).
# It then searches for SNPs in these regions from a file formatted:
#
#ZMLCHR10P1 10 3367 0.6916309672705971
#ZMLCHR10P2 10 3372 0.23592115510486514
#ZMLCHR10P3 10 20354 0.3317208805970703
#
# and exports all of the SNPs that lie in the designated regions in a .txt file.  


def main():
    file_name = input("Enter name of file containing regions to export: ") # region file (described above)
    file_name_w = input("Enter name of output file: ") # output file for snp list
    print() # blank line for easier reading

    file = open(file_name, "r")

    low = [] # empty lists for processing
    high = []
    snp = []
    pos = []
    pr = []

    ##### FIND REGIONS TO EXPORT FROM FILE #####
    
    for line in file: # file containing regions

        y = line

        check = y.split(',') # split each line at commas

        if check[0] == "Chr10": # which chromosome to use- CHANGE THIS
            low.append(str(check[1])) # add lower limit to list
            high.append(str(check[2])) # add upper limit to list

        else:
            pass # do nothing

    file.close() # close input file

    ##### GET SNP NAMES FROM SNP FILE #####
    
    file_name2 = input("Enter the name of file containing SNP names: ")
    file2 = open(file_name2, "r")

    for line in file2: # file containing SNP names for particular chromosome
        z = line

        p = z.split() # split at spaces 
        
        snp.append(str(p[0])) # append SNP name to list snp
        pos.append(str(p[2])) # append SNP position to list pos
        
    file2.close() 

    ##### FIND SNPS WITHIN GIVEN REGIONS #####

    for i in range(0,len(low)): # do this for the number of regions we need to find
        bound_l = eval(low.pop()) # assigns value of last entry in list to variable
        bound_h = eval(high.pop())# etc.

        for i in range(len(pos)): # do this for every SNP from SNP file

            p = eval(pos[i]) # change value at position [i] to a number for comparison in next line
            
            if bound_h >= p >= bound_l: # if the pos of the SNP is within the region of interest
                pr.append(str(snp[i])) # add to list to print later

            else:
                pass
    ##### WRITE RESULTS TO OUTPUT FILE #####
            
    out = open(file_name_w, "w") # open file to be written (change file in " ")
    out.write(' \n'.join(pr)) # write the SNP list to new file
    out.close()
                
            

    

main()
