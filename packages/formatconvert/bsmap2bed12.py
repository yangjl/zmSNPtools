#!/usr/bin/env python
from __future__ import division
from time import sleep
#import pybedtools
import sys, os
import argparse
import textwrap
import timeit

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
    
def read_write_BS_bed(infile_name, outfile_name):
    #temp = []
    with open(infile_name) as foo:
        totlines = len(foo.readlines())
    with open(infile_name, 'r') as infile, open(outfile_name, 'w') as outfile: 
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
            #temp.append([tokens[0], int(tokens[1]) - 1, tokens[1], tokens[3],tokens[4], tokens[2], 
            #             tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10] ])
                         
            outfile.write("\t".join(map(str, ([tokens[0], int(tokens[1]) - 1, tokens[1], tokens[3],tokens[4], 
            tokens[2], tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10] ]))) + "\n")
            i = i + 1
            sys.stdout.write('\r')
            # the exact output you're looking for:
            sys.stdout.write("Reading [ %d%% ] of the total [ %s ]" % (round(i/totlines*100, 0), totlines))
            sys.stdout.flush()

#get_loci_info(y[4:31])
def version():
    ver0 = """
    ##########################################################################################
    BSMAP to BED12
    Author: Jinliang Yang to my little girl Olivia
    purpose: convert BSMAP format to BED12 format
    --------------------------------
    
    updated: 6/1/2015
        o version 0.9
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
    parser.add_argument('-o', '--output', help='BED12', type=str)
    return parser
    #parser = get_parser()
    #parser.print_help()

def main():
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['input'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### cal running time ######
    st = timeit.default_timer()
   
    read_write_BS_bed(infile_name=args['input'], outfile_name=args['output'])

    et = timeit.default_timer()
    print("\n")
    print(">>> [ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print(">>> Job finished!")


if __name__ == '__main__':
    main()    


