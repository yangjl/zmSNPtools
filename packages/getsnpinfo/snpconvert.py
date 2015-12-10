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
    snpconvert version 0.1
    Jinliang Yang

    --------------------------------
    Convert selected SNP in Hapmap format into 0/1/2 format

    USAGE:
    ##########################################################################################
    """
    return v0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

def read_snpinfo(infofile):
    with open(infofile, 'r') as infile:
        header = infile.readline()
        header = header.split()

        if not "rs#" in header:
            warning("You must have 'rs#' column in your header!")
        if not "major" in header:
            warning("You must have 'major' column in your header!")
        if not "minor" in header:
            warning("You must have 'minor' column in your header!")

        snpinfo = {}
        for line in infile:
            tokens = line.split()
            for idx, item in enumerate(header):
                #setdefault function says "Get the value with this key, or
                ## if that key isn't there, add this value and then return it."
                snpinfo.setdefault(item,[]).append(tokens[idx])
    return snpinfo

##########################################################################################
def readfile_and_process(infile_name, outfile_name):
    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                pass
            elif line.startswith("rs#"):
                tokens = line.split()
                outfile.write("\t".join(["snpid", "major", "minor"]) + "\t" +
                              "\t".join(tokens[(args['start']-1):]) + "\n")
            else:
                tokens = line.split()
                if not tokens[0] in snpinfo['rs#']:
                    pass
                else:
                    idx = snpinfo['rs#'].index(tokens[0])
                    major = snpinfo['major'][idx]
                    minor = snpinfo['minor'][idx]

                    snptokens = tokens[(args['start']-1):]
                    newsnp = snp_recode(snptokens, major, minor)
                    ### write the first serveral columns:
                    outfile.write("\t".join([tokens[0], major, minor]) + "\t" + "\t".join(newsnp) + "\n")


def snp_recode(snptokens, major, minor):

    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
    if len(snpset) > 3:
        warning("Observed more than three alleles!", snpset)

    x = []
    for item in snptokens:
        if item == "N":
            x.append('3')
        elif item == major:
            x.append('0')
        elif item == minor:
            x.append('2')
        else:
            x.append('1')
    return x


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
    parser.add_argument('-a','--hmp', help='Hapmap SNP file', type=str)
    parser.add_argument('-i','--info', help='SNP information file', type=str)
    parser.add_argument('-s','--start', help='start cols (1-based) of the genotype', type=int)
    #parser.add_argument('-e','--end', help='end cols (1-based) of the genotype', type=int)
    #parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="N")
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)

    return parser
    #parser = get_parser()
    #parser.print_help()

if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['info'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### read in the input file ######
    st = timeit.default_timer()

    snpinfo = read_snpinfo(args['info'])
    print(">>> Loaded [",  "%.f" % len(snpinfo['rs#']), " ] SNP from SNP info file.")

    print(">>> Recoding and writing SNP ...")
    readfile_and_process(args['hmp'], args['output'])

    et = timeit.default_timer()

    print(">>> [ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print(">>> Job finished!")
