#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os

def version():
    v0 = """
    ##########################################################################################
    snpfrq version 3.0
    Jinliang Yang
    updated: August.8th.2017, for eMaize data
    --------------------------------

    new function: add maf and missingness filtration!
    compute SNP frq and loci missing rate from DSF
    ##########################################################################################
    """
    return v0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def readfile_and_process(infile_name, outfile_name):

    with open(infile_name, 'r') as infile:
    
        line1 = infile.readline()
        line1array = line1.split()
        if(line1array[0] == "snp"):
            wtype = 1
        elif(line1array[0] == "chr" and line1array[1] == "pos"):
            wtype = 2

        line2 = infile.readline()
        line2array = line2.split()
        if(len(str(line2array[start])) != 1):
            warning("SNP should be coded with one number, like:'0/1/2'!")

    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile:

        line1 = infile.readline()
        line1data = line1.split()

        outfile.write("\t".join(["snpid","sum", "c0", "c1", "c2"]) + "\n")
        
        for line in infile:
            tokens = line.split()
            #tokens = ["N" if x== args['missingcode'] else x for x in tokens]
            out = get_loci_info(tokens)

            ### print out the results
            if out:
                outfile.write("\t".join([ tokens[0], str(out['sum']), str(out['c0']), str(out['c1']), str(out['c2']) ]) + "\n")
                



def get_loci_info(tokens):

    snptokens = tokens[start:]
    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = set0
    
    info = {}
    
    if len(snpset) > 3:
        warning("SNPs should only be '0/1/2'!")
    else:
        info['sum'] = sum(int(i) for i in snptokens)
        info['c0'] = snptokens.count("0")
        info['c1'] = snptokens.count("1")
        info['c2'] = snptokens.count("2")

    return info

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
    parser.add_argument('-i','--input', help='input file', type=str)
    parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="N")
    
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)

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

    ##### read in the input file ######
    st = timeit.default_timer()
    print("Reading and writing SNP info ...")

    start = 1
    readfile_and_process(args['input'], args['output'])

    et = timeit.default_timer()

    print("[ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print("Job finished!")
