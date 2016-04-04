#!/usr/bin/env python
from __future__ import division
import re
import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO

def get_ctxnum(reffile):
    """
    Get the number of CG/CHG/CHH from a reference genome FASTA file
    """
    with open(reffile) as infile:
        fasta = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
        for chr in fasta:
            fasta[chr] = str(fasta[chr].seq).upper()
    num_cg = 0
    num_chg = 0
    num_chh = 0
    for chr in fasta:
        num_cg += len([match.start() for match in re.finditer(r'(?=(CG))', fasta[chr])])
        num_cg += len([match.start()-1 for match in re.finditer(r'(?<=(CG))', fasta[chr])])
        num_chg += len([match.start() for match in re.finditer(r'(?=(C[ACT]G))', fasta[chr])])
        num_chg += len([match.start()-1 for match in re.finditer(r'(?<=(C[AGT]G))', fasta[chr])])
        num_chh += len([match.start() for match in re.finditer(r'(?=(C[ACT][ACT]))', fasta[chr])])
        num_chh += len([match.start()-1 for match in re.finditer(r'(?<=([AGT][AGT]G))', fasta[chr])])
    return num_cg, num_chg, num_chh
            
def read_ctxcov(cgmapfile):
    """
    Read the column of coverage from CGmap file
    """
    cgmap = pd.read_table(cgmapfile, header=True, usecols=[3, 4, 5], names=['ctx', "ratio", 'cov'])
    cov_cg = cgmap[cgmap['ctx'] == 'CG']['cov']
    cov_chg = cgmap[cgmap['ctx'] == 'CHG']['cov']
    cov_chh = cgmap[cgmap['ctx'] == 'CHH']['cov']
    ratio_cg = cgmap[cgmap['ctx'] == 'CG']['ratio'].mean()
    ratio_chg = cgmap[cgmap['ctx'] == 'CHG']['ratio'].mean()
    ratio_chh = cgmap[cgmap['ctx'] == 'CHH']['ratio'].mean()
    tot_cg = len(cov_cg)
    tot_chg = len(cov_chg)
    tot_chh = len(cov_chh)
    cov_cg = cov_cg.mean()
    cov_chg = cov_chg.mean()
    cov_chh = cov_chh.mean()

    return [cov_cg, cov_chg, cov_chh,
            ratio_cg, ratio_chg, ratio_chh,
            tot_cg, tot_chg, tot_chh]


#get_loci_info(y[4:31])
def version():
    ver0 = """
    ##########################################################################################
    BSMAP to VCF v0.3
    Author: Jinliang Yang to my little girl Olivia
    purpose: convert BSMAP format to VCF format
    --------------------------------
    
    updated: 2/23/2016
        o VCFv4.2
        o A/T alleles and all PASS
        o tags => GT:DP:CC:CT:GL
    ##########################################################################################
    """
    return ver0
    
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
    parser.add_argument('-i','--input', help='MSMAP output', type=str)
    parser.add_argument('-o', '--output', help='Convert to VCFv4.2 format', type=str)
    
    parser.add_argument('-l', '--lower', help='default=0.3, lower value of the ratio', type=float, default=0.3)
    parser.add_argument('-u', '--upper', help='default=0.7, upper value of the ratio', type=float, default=0.7)
    parser.add_argument('-v', '--verbose', help='default=0, no progress', type=int, default=0)
    return parser
    #parser = get_parser()
    #parser.print_help()

def main():
    parser = get_parser()
    args = vars(parser.parse_args())
    
    if args['path'] is not None:
        os.chdir(args['path'])

    if args['input'] is not None:
        print(version())
        
        ##### cal running time ######
        st = timeit.default_timer()
        id = args['input'].split("_")
        
        write_VCF_meta(outfile_name=args['output'], sampleid=id[0])

        read_write_BS_VCF(infile_name=args['input'], outfile_name=args['output'], verbose=args['verbose'], lower=args['lower'], upper=args['upper'])

        et = timeit.default_timer()
        print("\n")
        print(">>> [ %.0f ] minutes of run time!" % ((et - st)/60))
        print(">>> Job finished!")
    else:
        print(parser.print_help())
        sys.exit(1)
    

if __name__ == '__main__':
    main()    
