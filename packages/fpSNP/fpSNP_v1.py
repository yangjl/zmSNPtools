#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os

def version():
    ver0 = """
    ##########################################################################################
    fpSNP version 1.0
    Author: Jinliang Yang
    purpose: find a set of population specific fingerprint SNPs
    --------------------------------
    
    updated: 4/20/2015, run for SeeDs and Ames inbred lines with GBS data
        o support multiple missing codes
        o snpid, chr, pos, ...
    updated: 1/28/2015, run for US GBS data
    updated: 10/2/2014, first piece of the code
    ##########################################################################################
    """
    return ver0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def checkFile(infile_snp):

    with open(infile_snp, 'r') as infile:
        ### check the first line
        line1 = infile.readline()
        line1a = line1.split()
        if(line1a[0] != "snpid" or line1a[1] != "chr" or line1a[2] != "pos"):
            warning("!!! DSF3+: ['snpid', 'chr', 'pos'] should be the first 3 cols in the header")
        else:
            print("--- input genotype file header OK!")

        ### check the 2nd line
        line2 = infile.readline()
        line2a = line2.split()
        for gi in line2a[3:]:
            if(not gi.isdigit()):
                warning("!!! group coding should be integer")
        else:
            print("--- input group coding OK!")

        ### check the 3rd line
        line3 = infile.readline()
        line3a = line3.split()
        for snpi in line3a[3:]:
            if(len(str(snpi)) != 1):
                warning("!!! SNP should be coded with one latter, like:'A T C G or - +'!")
        else:
            print("--- input genotype file coding format OK!")

#f='/Users/yangjl/Desktop/simusnps.txt'
#checkFile(f)

##########################################################################################
def readfile_and_process(infile_snp):

    with open(infile_snp, 'r') as infile:

        line1 = infile.readline()
        ### idx the groups
        line2 = infile.readline()
        line2a = line2.split()

        set0 = set(line2a[3:])
        if 'N' in set0:
            set0.remove('N')
        groups = list(set0)
        #print(groups)
        if(len(groups) !=2):
            warning("Oops, you have more than two groups")
        else:
            idx1 = [i for i, j in enumerate(line2a) if j == groups[0]]
            #print(idx1)
            idx2 = [i for i, j in enumerate(line2a) if j == groups[1]]
            #print(idx2)

        for line in infile:
            tokens = line.split()
            ### change multiple missing codes to N
            mcode = list(args['missingcode'])
            for amcode in mcode:
                tokens = ["N" if x== amcode else x for x in tokens]
            out = get_loci_info(tokens, idx= range(3,len(tokens)), MA="N")
            
            if 'minor' in out:
                out1 = get_loci_info(tokens, idx = idx1, MA = out['minor'])
                out2 = get_loci_info(tokens, idx = idx2, MA = out['minor'])
                snpmaf.append([tokens[0], tokens[1], tokens[2], out['major'],out['minor'], out['maf'], out['missing'],
                           out1['maf'],out1['missing'], out2['maf'], out2['missing']])
#readfile_and_process(f)

def write_maf():
    ### print out the results
    with open(args['output'], 'w') as outfile:
        outfile.write("\t".join(['snpid', "chr", "pos", "major", "minor", "maf0", "missing0", "maf1", "missing1", "maf2", "missing2"]) + "\n")
        for mysnp in snpmaf:
            outfile.write("\t".join(str(v) for v in mysnp) + "\n")
    print("--- MAF file wrote to [ {} ]".format(args['output']))

def write_prob_snp():
    """
    Write the problem snp to the input.prob file:
    """
    output = args['output'].split('.')[0]
    output1 = ".".join([output, "one"])
    outfile1 = open(output1, 'w')
    for prob_snp1 in prob1:
        outfile1.write('\t'.join(prob_snp1) + '\n')
    outfile1.close()
    print("--- loci with no polymorphism wrote to [ {} ]".format(output1))

    output3 = ".".join([output, "multiple"])
    outfile3 = open(output3, 'w')
    for prob_snp3 in prob3:
        outfile3.write('\t'.join(prob_snp3) + '\n')
    outfile3.close()
    print("--- loci with multiple alleles wrote to [ {} ]".format(output3))



### function to calculate maf for a given set of plants
### if Minor allele (MA) = "N", minor allele will be determined by the program
def get_loci_info(tokens, idx=[], MA="N"):

    snptokens = [tokens[i] for i in idx]
    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    info = {}
    if len(snpset) == 1:
        prob1.append(tokens)
        if snpset[0] == MA:
            info['maf'] = 1
        else:
            info['maf'] = 0

        c1 = snptokens.count("N")
        info['missing'] = round(c1/len(snptokens), 5)

    if len(snpset) > 2:
        prob3.append(tokens)
    if len(snpset) == 2:
        if MA == "N":
            c1 = snptokens.count(snpset[0])
            c2 = snptokens.count(snpset[1])
            info['missing'] = round((len(snptokens) - c1 - c2)/len(snptokens), 5)

            if c1 >= c2:
                info['major'] = snpset[0]
                info['minor'] = snpset[1]
                info['maf'] = round(c2/(c1+c2), 5)
            else:
                info['major'] = snpset[1]
                info['minor'] = snpset[0]
                info['maf'] = round(c1/(c1+c2), 5)
        elif MA != "N":
            c1 = snptokens.count(MA)
            cn = snptokens.count("N")
            info['maf'] = round(c1/(len(snptokens) - cn), 5)
            info['missing'] = round(cn/len(snptokens), 5)
    return info
####
##########################################################################################
#get_loci_info(y[4:31])
    
def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(version())
        )

    # positional arguments:
    #parser.add_argument('query', metavar='QUERY', type=str, nargs='*', \
    #                    help='the question to answer')

    # optional arguments:
    parser.add_argument('-p', '--path', help='the path of the input files', \
                        nargs='?', default=os.getcwd())
    parser.add_argument('-i','--input', help='input file: 3+ Density SNP Format (DSF), [snpid, chr, pos, plant1, ...]', type=str)

    parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="-0")

    parser.add_argument('-o', '--output', help='output files, such as chr1_merged', type=str)

    return parser
    #parser = get_parser()
    #parser.print_help()

if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['input'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### cal running time ######
    st = timeit.default_timer()
    prob1=[]
    prob3=[]
    snpmaf=[]
    print(">>> Reading and writing SNP info ...")
    ### checking the input file
    checkFile(args['input'])
    ### read the input file
    readfile_and_process(args['input'])

    #### writing the output file
    print(">>> Writing SNPs with multiple alleles ...")
    write_prob_snp()
    write_maf()

    et = timeit.default_timer()

    print(">>> [ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print(">>> Job finished!")
