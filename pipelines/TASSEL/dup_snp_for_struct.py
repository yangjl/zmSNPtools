# Duplicate SNP Data for STRUCTURE v.2
# Chris Fiscus
# 2/3/14

# This script copies SNP data after it has been transposed so it can be used in structure.
# If file already has two lines for every individual then no need to use this script. 

def main():
    f_r = open("reformatted.txt", "r") # file to read from
    f_w = open("outfile.txt", "w") # file to write to

    x = 0 # used to skip first line in file

    for line in f_r: # do this for every line in the file
        
        if x == 0: # first line in file
            line = line.replace("rs", "") # remove string "rs" 
            f_w.write(line) # copy line from input file to output file

        else: # every line after first line
            f_w.write(line) # write lines twice (for structure format)
            f_w.write(line)

        x = x + 1 # line number 



    f_r.close() # close files so they can be used
    f_w.close()

main()
