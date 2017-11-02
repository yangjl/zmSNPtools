#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os
import pandas as pd

def version():
    v0 = """
    ##########################################################################################
    snpfrq2 version 4.0
    Jinliang Yang
    updated: Octo.24.2017, count SNP frq for g1 and g2
    --------------------------------
    
    This version of snpfrq can determine the major/minor allele, calculate frq for two subgroups.
    new function: add maf and missingness filtration!
    compute SNP frq and loci missing rate from DSF
    ##########################################################################################
    """
    return v0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()


def get_groups(groupfile):
    df = pd.read_csv(groupfile)
    ug = df.group.unique()
    if(len(ug) != 2):
        warings("group file should have header: group (two levels) and pid.")
    g1 = df[df.group == ug[0]].pid
    g2 = df[df.group == ug[1]].pid
    gout = {'g1':list(map(str, g1)), 'g2':list(map(str, g2))}
    return gout

##########################################################################################
def readfile_and_process(infile_name, outfile_name, gout):

    with open(infile_name, 'r') as infile:
    
        line1 = infile.readline()
        line1array = line1.split()
        
        # determine the idx of the two groups
        idx1 = []
        for i in gout['g1']:
            idx1.append(line1array.index(i))
        idx2 = []
        for j in gout['g2']:
            idx2.append(line1array.index(j))

        line2 = infile.readline()
        line2array = line2.split()
        if(len(str(line2array[start])) != 1):
            warning("SNP should be coded with one latter, like:'A T C G or - +'!")

    myout = args['output'].split('.')[0]
    myout = ".".join([myout, "frq"])
    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile:

        line1 = infile.readline()
        line1data = line1.split()

        outfile.write("\t".join(["snpid","major","minor","MAF","missing", "majorg1", "minorg1", "majorg2", "minorg2"])  + "\n")
        
        for line in infile:
            tokens = line.split()
            tokens = ["N" if x == args['missingcode'] else x for x in tokens]
            tokens = ["N" if x not in ('A', 'T', 'C', 'G') else x for x in tokens]
            #idx1 and idx2 would be the idx of the two groups.
            out = get_loci_info(tokens, idx1, idx2)

            ### print out the results
            if out:
                outfile.write("\t".join([out['snpid'], out['major'], out['minor'], str(out['maf']), str(out['missing']),
                str(out['cmajor1']),str(out['cminor1']), str(out['cmajor2']), str(out['cminor2'])] )  + "\n")   
                

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

    output3 = ".".join([output, "mul"])
    outfile3 = open(output3, 'w')
    for prob_snp3 in prob3:
        outfile3.write('\t'.join(prob_snp3) + '\n')
    outfile3.close()

def get_loci_info(tokens, idx1, idx2):

    snptokens = tokens[start:end]
    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    info = {}
    if len(snpset) == 1:
        prob1.append(tokens)
    if len(snpset) > 2:
        prob3.append(tokens)
    elif len(snpset) == 2:
        c1 = snptokens.count(snpset[0])
        c2 = snptokens.count(snpset[1])
   
        if c1 >= c2:
            info['major'] = snpset[0]
            info['minor'] = snpset[1]
            info['maf'] = round(c2/(c1+c2),3)
            info['cmajor1'] = [tokens[i] for i in idx1].count(snpset[0]) #[L[i] for i in Idx]
            info['cminor1'] = [tokens[i] for i in idx1].count(snpset[1])
            info['cmajor2'] = [tokens[i] for i in idx2].count(snpset[0])
            info['cminor2'] = [tokens[i] for i in idx2].count(snpset[1])
        else:
            info['major'] = snpset[1]
            info['minor'] = snpset[0]
            info['maf'] = round(c1/(c1+c2),3)
            info['cmajor1'] = [tokens[i] for i in idx1].count(snpset[1])
            info['cminor1'] = [tokens[i] for i in idx1].count(snpset[0])
            info['cmajor2'] = [tokens[i] for i in idx2].count(snpset[1])
            info['cminor2'] = [tokens[i] for i in idx2].count(snpset[0])
        info['missing'] = round((len(snptokens) - c1 - c2)/len(snptokens),3)
        info['snpid'] = tokens[0]

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
    parser.add_argument('-g','--group', help='a file with group information', type=str)
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)
    
    parser.add_argument('-s','--start', help='start cols (1-based) of the genotype', type=int)
    parser.add_argument('-e','--end', help='end cols (1-based) of the genotype', type=int)
    parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="N")

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
    prob1=[]
    prob3=[]
    print("Reading and writing SNP info ...")

    start = args['start'] -1
    end = args['end']
    gout = get_groups(groupfile=args['group'])
    readfile_and_process(args['input'], args['output'], gout)

    #print("Writing SNPs with multiple alleles ...")
    #write_prob_snp()

    et = timeit.default_timer()

    print("[ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print("Job finished!")
