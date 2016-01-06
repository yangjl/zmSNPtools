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
    getvcf version 0.1
    Jinliang Yang
    Jan 6th, 2015, extract data from VCF

    --------------------------------
    USAGE: getvcf -i test.vcf -n snpid.txt -o out3.txt
    --------------------------------
    Note:
    ##########################################################################################
    """
    return v0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def read_snpid(snpfile_name):
    with open(snpfile_name, 'r') as snpfile:
        line1 = snpfile.readline()

        snpid = []
        for line in snpfile:
            tokens = line.split()
            snpid.append(tokens[0])
    return snpid

def readfile_and_process(infile_name, outfile_name, snpid):

    nsnp = 0
    p = 1
    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile:
        for line in infile:
            ### print progress
            if nsnp/args['barlen'] > p:
                p = p +1
                print("###>>> processed [ ", "%s" % nsnp, " ]  SNPs")

            if line.startswith("##"):
                pass
            elif line.startswith("#CHROM"):
                hd = line.split()
                outfile.write("CHROM" + "\t")
                outfile.write("\t".join(hd[1:]) + "\n")
            else:
                tokens = line.split()

                if tokens[2] in snpid:
                    nsnp = nsnp + 1
                    ### write the first serveral columns:
                    outfile.write("\t".join(tokens[0:(args['start']-1)]) )
                    for info in tokens[(args['start']-1):]:
                        DP = info.split(":")
                        outfile.write("\t" + DP[args['element']])
                    outfile.write("\n")


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(version())
        )
    # optional arguments:
    parser.add_argument('-p', '--path', help='the path of the input files', \
                        nargs='?', default=os.getcwd())
    parser.add_argument('-i','--input', help='input file', type=str)
    parser.add_argument('-n','--snpid', help='txt file of snpid, chr, pos', type=str)
    parser.add_argument('-s','--start', help='start cols (1-based) of the genotype, default=10', type=int, default=10)
    parser.add_argument('-e','--element', help='element to extract from the slot, default=2', type=int, default=2)
    parser.add_argument('-b','--barlen', help='number to print, default=10000', type=int, default=10000)
    #parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="N")
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['input'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### read in the input file ######
    st = timeit.default_timer()
    snpid = read_snpid(snpfile_name = args['snpid'])
    print("###>>> Loaded [ %s ] SNPs!" % len(snpid))

    ##### start to process data ######
    readfile_and_process(infile_name=args['input'], outfile_name=args['output'], snpid=snpid)

    et = timeit.default_timer()

    print("###>>> Run time [ %.0f ] minutes!" % ((et - st)/60))
    print("###>>> Job finished!")
