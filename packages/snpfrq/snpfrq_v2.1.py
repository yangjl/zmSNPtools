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
    snpfrq version 2.1
    Jinliang Yang
    updated: April 8th, 2015, for SAM SNPs
    --------------------------------
    compute SNP frq and loci missing rate from 'BED+' format
    USAGE: ./snpfrq_v2.1.py -i test_hmp1_chr1.bed -s 5 -m "-N+" -a 0 -b 1 -c 2 -o test.out

    --------------------------------
    new feature: 1. support multiple missing codes
    new feature: 2. add maf and missingness filtration!

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
        if(line1array[3] == "snpid"):
            wtype = 1
        elif(line1array[0] == "chr" and line1array[1] == "pos"):
            wtype = 2

        line2 = infile.readline()
        line2array = line2.split()
        if(len(str(line2array[start])) != 1 and len(str(line2array[start])) != 3):
            warning("SNP should be coded with one letter, like:'A T C G or - +'!" or 'A A, T T, ...')

    myout = args['output'].split('.')[0]
    myout = ".".join([myout, "flt"])
    with open(infile_name, 'r') as infile, open(outfile_name, "w") as outfile, open(myout, "w") as outfile2:

        line1 = infile.readline()
        line1data = line1.split()
        
        
        if args['outmode'] == 1:
            outfile.write("\t".join(["snpid","major","minor","MAF","missing"]) + "\t" + "\t".join(line1data[start:])  + "\n")
            outfile2.write("\t".join(["snpid","major","minor","MAF","missing"]) + "\t" + "\t".join(line1data[start:]) + "\n")
        elif args['outmode'] == 2:
            outfile.write("\t".join(["snpid","major","minor","MAF","missing"]) + "\n")
            outfile2.write("\t".join(["snpid","major","minor","MAF","missing"]) + "\n")


        for line in infile:
            tokens = line.split()

            ### change the missing code to N
            mcode = list(args['missingcode'])
            for amcode in mcode:
                tokens = ["N" if x== amcode else x for x in tokens]

            ### get information for each locus
            out = get_loci_info(tokens)

            ### print out the results
            if out:
                if wtype == 1 and args['outmode'] == 1:
                    if out['maf'] >= args['maf'] and out['missing'] <= args['mr']:
                        outfile.write("\t".join([ tokens[3], out['major'], out['minor'], str(out['maf']), str(out['missing']) ]) + \
                          "\t" + "\t".join(tokens[start:]) + "\n")
                    else:
                        outfile2.write("\t".join([tokens[3], out['major'], out['minor'], str(out['maf']), str(out['missing'])]) + "\t" \
                           + "\t".join(tokens[start:]) + "\n")
                if wtype == 1 and args['outmode'] == 2:
                    if out['maf'] >= args['maf'] and out['missing'] <= args['mr']:
                        outfile.write("\t".join([ tokens[3], out['major'], out['minor'], str(out['maf']), \
                        str(out['missing']) ]) + "\n")
                    else:
                        outfile2.write("\t".join([tokens[3], out['major'], out['minor'], str(out['maf']), \
                        str(out['missing'])]) + "\n")           
                elif wtype == 2:
                    if out['maf'] >= args['maf'] and out['missing'] <= args['mr']:
                        outfile.write("_".join([ tokens[0],tokens[1] ]) + "\t" + \
                          "\t".join([out['major'], out['minor'], str(out['maf']), str(out['missing'])]) + \
                          "\t" + "\t".join(tokens[start:]) + "\n")
                    else:
                        outfile2.write("_".join([ tokens[0],tokens[1] ]) + "\t" + \
                           "\t".join([out['major'], out['minor'], str(out['maf']), str(out['missing'])]) + \
                           "\t" + "\t".join(tokens[start:]) + "\n")


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

def get_loci_info(tokens):

    snptokens = tokens[start:]
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
        else:
            info['major'] = snpset[1]
            info['minor'] = snpset[0]
            info['maf'] = round(c1/(c1+c2),3)
        info['missing'] = round((len(snptokens) - c1 - c2)/len(snptokens),3)

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
    parser.add_argument('-s','--start', help='start cols (1-based) of the genotype', type=int)
    #parser.add_argument('-e','--end', help='end cols (1-based) of the genotype', type=int)
    parser.add_argument('-m','--missingcode', help='code for missingness', type=str, default="N")
    parser.add_argument('-a','--maf', help='cutoff of minor allele freq', type=float, default=0)
    parser.add_argument('-b','--mr', help='cutoff of missing rate', type=float, default=1)
    parser.add_argument('-c','--outmode', help='mode of output', type=float, default=2)
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
    prob1=[]
    prob3=[]
    print("Reading and writing SNP info ...")

    start = args['start'] -1
    #end = args['end']
    readfile_and_process(args['input'], args['output'])

    print("Writing SNPs with multiple alleles ...")
    write_prob_snp()

    et = timeit.default_timer()

    print("[ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print("Job finished!")
