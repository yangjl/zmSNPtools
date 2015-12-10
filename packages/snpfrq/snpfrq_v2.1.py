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
    snpfrq version 3.2
    Jinliang Yang
    updated: Dec 8th, 2015, for HapMap format GBS SNPs
    updated: add IUPAC code
    --------------------------------
    Compute SNP frq and loci missing rate from 'HapMap or BED+' format

    USAGE: snpfrq -i test_hmp1_chr1.dsf -s 5 -m "N" -o out.txt
    --------------------------------
    Note:  maf=0 for non-variant sites. maf=1 for sites with multiple alleles.
    missing sets to -9 for these sites.
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
            elif line.startswith("rs#"):
                line1array = line.split()
                outfile.write("\t".join(line1array[0:(args['start']-1)]) + "\t" + "\t".join(["major", "minor", "MAF", "missing"]) + "\n")
            else:
                tokens = line.split()
                ### write the first serveral columns:
                outfile.write("\t".join(tokens[0:(args['start']-1)]) + "\t")
                snptokens = tokens[(args['start']-1):]
                alleles = tokens[1]

                ### change the missing code to N
                #mcode = list(args['missingcode'])
                #for amcode in mcode:
                #    snptokens = ["N" if x== amcode else x for x in snptokens]

                ### replace with IUPAC
                snptokens = IUPAC(snptokens, alleles)

                ### get information for each locus
                out = get_loci_info(snptokens)
                outfile.write("\t".join([ out['major'], out['minor'], str(out['maf']), str(out['missing']) ]) + "\n")


def get_loci_info(snptokens):

    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    info = {}
    if len(snpset) < 2:
        info['major'] = ''.join(snpset)
        info['minor'] = ''.join(snpset)
        info['maf'] = 0
        info['missing'] = -9
    elif len(snpset) > 2:
        info['major'] = ''.join(snpset)
        info['minor'] = ''.join(snpset)
        info['maf'] = 1
        info['missing'] = -9
    elif len(snpset) == 2:
        c1 = snptokens.count(snpset[0])
        c2 = snptokens.count(snpset[1])
   
        if c1 >= c2:
            info['major'] = snpset[0]
            info['minor'] = snpset[1]
            info['maf'] = round(c2/(c1+c2), 5)
        else:
            info['major'] = snpset[1]
            info['minor'] = snpset[0]
            info['maf'] = round(c1/(c1+c2), 5)
        info['missing'] = round((len(snptokens) - c1 - c2)/len(snptokens), 5)

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
    #parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="N")
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
