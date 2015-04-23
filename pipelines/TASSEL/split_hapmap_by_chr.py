# Format Hapmap
# Chris Fiscus
# 2/13/14
# A script that sorts a hapmap file by chromosome and prints each chromosome into a separate file 

# Before running this script, you need to change extension of .hapmap file to .txt


def main():

    file = open("SNP_rand_10k.txt", "r") # open hapmap file to be processed (change file in " ")

    x1 = [] # chromosome 1 
    x2 = [] # chromosome 2
    x3 = [] # chromosome 3 
    x4 = [] # etc. 
    x5 = []
    x6 = []
    x7 = []
    x8 = []
    x9 = []
    x10 = []
    xp = [] # header of input file 

    
    for line in file: # do this for every line in the file

        y = line

        check = y.split() # split line at spaces

        if check[2] == "chr1":   # sort lines by chromosome 
            x1.append(str(line)) # write line to list 

        elif check[2] == "chr2": # etc. 
            x2.append(str(line))

        elif check[2] == "chr3":
            x3.append(str(line))

        elif check[2] == "chr4":
            x4.append(str(line))

        elif check[2] == "chr5":
            x5.append(str(line))

        elif check[2] == "chr6":
            x6.append(str(line))

        elif check[2] == "chr7":
            x7.append(str(line))

        elif check[2] == "chr8":
            x8.append(str(line))

        elif check[2] == "chr9":
            x9.append(str(line))

        elif check[2] == "chr10":
            x10.append(str(line))

        else:  # header of file 
            xp.append(str(line))

    file.close() # close input file

    out = open("header.txt","w") # open file to be written (for file header)
    out.write(''.join(xp)) # write the contents of list to new file
    out.close # close out file 

    out = open("chr1.txt", "w") # open file to be written (change file in " ")
    out.write(''.join(x1)) # write the list to new file
    x1[:] = [] # empty list
    out.close # close out file

    out = open("chr2.txt", "w") # etc. 
    out.write(''.join(x2))
    x2[:] = []
    out.close

    out = open("chr3.txt","w")
    out.write(''.join(x3))
    x3[:] = []
    out.close

    out = open("chr4.txt","w")
    out.write(''.join(x4))
    x4[:] = []
    out.close

    out = open("chr5.txt","w")
    out.write(''.join(x5))
    x5[:] = []
    out.close

    out = open("chr6.txt","w")
    out.write(''.join(x6))
    x6[:] = []
    out.close

    out = open("chr7.txt","w")
    out.write(''.join(x7))
    x7[:] = []
    out.close

    out = open("chr8.txt","w")
    out.write(''.join(x8))
    x8[:] = []
    out.close

    out = open("chr9.txt","w")
    out.write(''.join(x9))
    x9[:] = []
    out.close

    out = open("chr10.txt","w")
    out.write(''.join(x10))
    x10[:] = []
    out.close

 

main()
