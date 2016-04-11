#!/usr/bin/env python
from __future__ import division
from time import sleep, gmtime, strftime
#import pybedtools
import sys, os
import argparse
import textwrap
import timeit
import math

##########################################################################################
# http://stackoverflow.com/questions/9268652/progress-bar-for-reading-lines-in-text-file
class SimpleProgressBar(object):
    def __init__(self, maximum, state=0):
        self.max = maximum
        self.state = state

    def _carriage_return(self):
        sys.stdout.write('\r')
        sys.stdout.flush()

    def _display(self):
        stars = ''.join(['*'] * self.state + [' '] * (self.max - self.state))
        print '[{0}] {1}/{2}'.format(stars, self.state, self.max),
        self._carriage_return()

    def update(self, value=None):
        if not value is None:
            self.state = value
        self._display()

#spb = SimpleProgressBar(10)
#for i in range(0, 11):
#    time.sleep(.3)
#    spb.update(i)
def write_VCF_meta(outfile_name, sampleid):
    
    with open(outfile_name, 'w') as outfile: 
        ### check the first line
        outfile.write("##fileformat=VCFv4.2" + "\n")
        outfile.write("##fileDate=" + strftime("%m-%d-%Y", gmtime()) + "\n")
        outfile.write("##source=bsmap2VCF.py" + "\n")
        outfile.write("##reference=AGPv2" + "\n")
        
        outfile.write('##INFO=<ID=CO,Number=1,Type=String,Description="Context: CG, CHG and CHH">' + "\n")
        outfile.write('##INFO=<ID=ST,Number=1,Type=String,Description="Strand">' + "\n")

        outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, 1 methylated and 0 non-methylated">' + "\n")
        outfile.write('##FORMAT=<ID=RA,Number=1,Type=Float,Description="Methylation Ratio, column ratio">' + "\n")
        
        outfile.write('##FORMAT=<ID=CC,Number=1,Type=Integer,Description="C Count column C_count">' + "\n")
        outfile.write('##FORMAT=<ID=CT,Number=1,Type=Integer,Description="CT Count=Depth, column CT_count">' + "\n")
        outfile.write('##FORMAT=<ID=GL,Number=3,Type=Float,Description="Genotype likelihood for 00/01/11">' + "\n")
        
        outfile.write("\t".join(map(str, (["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sampleid ])))  + "\n")      

#######            
def read_write_BS_VCF(infile_name, outfile_name,  verbose, lower, upper):
    #temp = []
    if(verbose > 0):
        with open(infile_name) as foo:
            totlines = len(foo.readlines())
    with open(infile_name, 'r') as infile, open(outfile_name, 'a') as outfile: 
        ### check the first line
        line1 = infile.readline()
        line1a = line1.split()
        if(line1a[0] != "chr" or line1a[1] != "pos" or line1a[2] != "strand" or line1a[3] != "context"):
            warning("!!! input file: ['chr', 'pos', 'strand', 'context'] should be the first 4 cols in the header")
        else:
            print("--- input genotype file header OK!")
        
        #spb = SimpleProgressBar(totline)
        i = 0
        for line in infile:
            tokens = line.split()
            #1. chrom, 2.start, 3.end, 4.name=context, 5, score, 6. strand, 7-12
            ###>>>
            outfile.write(tokens[0] + "\t" + tokens[1] + "\t" + tokens[0] + "_" + tokens[1] + "\t" + "C" + "\t" + "T" + "\t")
            outfile.write(str(int(float(tokens[5]))) + "\t")
            #if(int(float(tokens[5])) < 5):
            #    outfile.write("q5" + "\t")
            #else:
            outfile.write("PASS" + "\t")
            ### INFO
            outfile.write("CO=" + tokens[3] + ";")
            outfile.write("ST=" + tokens[2] + "\t")
            
            ### FORMAT
            outfile.write("GT:RA:CC:CT:GL" + "\t")
            if(float(tokens[4]) < lower and float(tokens[4]) >= 0):
                outfile.write("0/0:" + tokens[4] + ":" + tokens[6] + ":" + tokens[7] + ":")
            elif(float(tokens[4]) > upper and float(tokens[4]) <= 1):
                outfile.write("1/1:" + tokens[4] + ":" + tokens[6] + ":" + tokens[7] + ":")
            else:
                outfile.write("0/1:" + tokens[4] + ":" + tokens[6] + ":" + tokens[7] + ":")
            
            q = float(tokens[6])/float(tokens[7])
            if(q == 0):
                q = 1/(10**5)
            elif(q == 1):
                q = 1 - 1/(10**5)
            p = 1 - q
            outfile.write(str(round(math.log10(p**2),2)) + ",")
            outfile.write(str(round(math.log10(2*p*q),2)) + ",")
            outfile.write(str(round(math.log10(q**2),2)) + "\n")
            
            #### print progress
            if(verbose > 0):
                i = i + 1
                sys.stdout.write('\r')
                # the exact output you're looking for:
                sys.stdout.write("Reading [ %d%% ] of the total [ %s ]" % (round(i/totlines*100, 0), totlines))
                sys.stdout.flush()


        
#get_loci_info(y[4:31])
def version():
    ver0 = """
    ##########################################################################################
    BSMAP to VCF v0.4
    Author: Jinliang Yang to my little girl Olivia
    purpose: convert BSMAP format to VCF format
    --------------------------------
    
    GT:RA:CC:CT:GL
    => GT: genotype, 0/0 unmethylated; 1/1 methylated; 0/1 partially methylated. Determined by l=0.3 and u=0.7.
    => RA: methylation ration, CC/CT, column ratio.
    => CC: C_count, methylated C.
    => CT: CT_count, total number of C and T.
    => GL: Genotype likelihood, p=ratio, q=1-p. Then computed using p^2, 2pq and q^2.
    
    updated: 2/23/2016
        o VCFv4.2
        o C/T alleles and all PASS
        o tags => GT:RA:CC:CT:GL
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


