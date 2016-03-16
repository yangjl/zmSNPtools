#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from collections import Counter
import sys
import argparse
import textwrap
import timeit
import os
import operator

def version():
    v0 = """
    ##########################################################################################
    gerpinfo version 0.1
    Jinliang Yang
    ##########################################################################################
    """
    return v0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def readfile_and_process(infile_name, outfile_name):
    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                pass
            elif line.startswith("V1"):
                line1array = line.split()
                outfile.write("\t".join(line1array[0:(args['start']-1)]) + "\t" +
                              "\t".join(["major", "minor", "majorc","minorc", "missing"]) + "\n")
            else:
                tokens = line.split()
                ### write the first serveral columns:
                outfile.write("\t".join(tokens[0:(args['start']-1)]) + "\t")
                snptokens = tokens[(args['start']-1):]
                #alleles = tokens[1]

                ### change the missing code to N
                mcode = list(args['missingcode'])
                for amcode in mcode:
                    snptokens = ["N" if x== amcode else x for x in snptokens]

                ### replace with IUPAC
                #snptokens = IUPAC(snptokens, alleles)

                ### get information for each locus
                out = get_loci_info(snptokens)
                outfile.write("\t".join([ out['major'], out['minor'], str(out['majorc']),str(out['minorc']), str(out['missing']) ]) + "\n")


def get_loci_info(snptokens):
    #snptokens = ['A', 'N', 'N', 'A', 'A', 'C']
    tot = len(snptokens)
    cN = 0
    c = Counter(snptokens)
    if 'N' in c:
        cN = c['N']
        c.pop('N')
    snpset = c

    info = {}
    info['missing'] = round(cN/tot, 5)
    if len(snpset) < 2:
        info['major'] = snpset.keys()[0]
        info['minor'] = 'N'
        info['majorc'] = snpset.values()[0]
        info['minorc'] = 0

    elif len(snpset) >= 2:
        info['major'] = snpset.keys()[0]
        info['minor'] = snpset.keys()[1]
        info['majorc'] = snpset.values()[0]
        info['minorc'] = snpset.values()[1]

    return info

### note this is hapmap version IUPAC, + actually= ATCG/-
def IUPAC(snptokens, alleles):

    alleles = alleles.split("/")
    if '-' in alleles:
        alleles.remove('-')

    x = []
    for index, item in enumerate(snptokens):
        if item == "N":
            x.extend(["N", "N"])
        elif item == "A":
            x.extend(["A", "A"])
        elif item == "C":
            x.extend(["C", "C"])
        elif item == "G":
            x.extend(["G", "G"])
        elif item == "T":
            x.extend(["T", "T"])
        elif item == "R":
            x.extend(["A", "G"])
        elif item == "Y":
            x.extend(["C", "T"])
        elif item == "S":
            x.extend(["C", "G"])
        elif item == "W":
            x.extend(["A", "T"])
        elif item == "K":
            x.extend(["G", "T"])
        elif item == "M":
            x.extend(["A", "C"])
        elif item == "+":
            x.extend(["+", "+"])
        elif item == "0":
            if len(alleles) == 1:
                x.extend([alleles[0], "-"])
            else:
                x.extend(["+", "-"])
        elif item == "-":
            x.extend(["-", "-"])
        else:
            warning("detected non-IUPAC genotype", item)
    return(x)

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
    parser.add_argument('-s','--start', help='start cols (1-based) of the genotype', type=int)
    #parser.add_argument('-e','--end', help='end cols (1-based) of the genotype', type=int)
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
    print(">>>Reading and writing SNP info ...")
    readfile_and_process(args['input'], args['output'])

    et = timeit.default_timer()

    print(">>>[ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print(">>>Job finished!")
