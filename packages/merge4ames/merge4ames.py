#!/usr/bin/env python

from collections import Counter
from __future__ import print_function
import sys
import argparse
import textwrap
import os

def version():
    v1 = """
    ##########################################################################################
    merge4ames version 1.0
    Jinliang Yang
    updated: 2.11.2014
    --------------------------------

    Merge replicated samples of the Ames panel GBS data to produce consensus SNP callings.
    When comparing two genotypes, only replace the genotype in the original list
    if it is missing and the second genotype is not or make it missing if the two genotypes
    do not agree!
    ##########################################################################################
    """
    return v1

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
    parser.add_argument('--hmp1', help='path of the input hmp1 data', type=str)
    parser.add_argument('--hmp2', help='path of the input hmp2 data', type=str)
    parser.add_argument('--rnaseq', help='path of the input rna-seq data', type=str)
    parser.add_argument('-f', '--infiles', help='a list of merging files', type=str)
    parser.add_argument('-o', '--output', help='output files, like chr1_merged', type=str)
    parser.add_argument('-m','--mode',
                        help='''[--mode=0 default], consensus, otherwise N;
                        1, keep hmp1;
                        2, keep hmp2;
                        3, keep RNA-seq;
                        4, vote the most freq call.
                        ''',
                        choices=[0,1,2,3,4], default=0, type=int)
    return parser

### use print_function
def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

##########################################################################################
def read_idx(ames_idx_file):
    temp = []
    infile = open(ames_idx_file, 'r')
    next(infile)
    for line in infile:
        line = line.strip()
        token = line.split('\t')
        temp.append(token)
    infile.close()
    return temp # return a list of list ['PI405705', '16 17 18 19 20 21', '6']

#The GBS genotype files in HapMap format (*.hmp.txt.gz) can be opened directly in TASSEL in
#their compressed (*.gz) state. To save space, we use single letters to encode unphased,
#diploid genotypes, where:
#A = A/A
#C = C/C
#G = G/G
#T = T/T
#K = G/T
#M = A/C
#R = A/G
#S = C/G
#W = A/T
#Y = C/T
#N = missing
def read_snp(snpfile):
    temsnp = []
    countN = 0
    countATCG = 0
    countY = 0
    infile = open(snpfile, 'r')
    
    fline = infile.readline()
    nms = fline.split()

    for line in infile:
        line = line.strip()
        token = line.split()
        for i in xrange(11, len(token)):
            if token[i] == 'A':
                token[i] = 'A A'
                countATCG += 1
            elif token[i] == 'C':
                token[i] = 'C C'
                countATCG += 1
            elif token[i] == 'G':
                token[i] = 'G G'
                countATCG += 1
            elif token[i] == 'T':
                token[i] = 'T T'
                countATCG += 1
            elif token[i] == 'K':
                token[i] == 'G T'
                countY += 1
            elif token[i] == 'M':
                token[i] = 'A C'
                countY += 1
            elif token[i] == 'R':
                token[i] = 'A G'
                countY += 1
            elif token[i] == 'S':
                token[i] = 'C G'
                countY += 1
            elif token[i] == 'W':
                token[i] = 'A T'
                countY += 1
            elif token[i] == 'Y':
                token[i] = 'C T'
                countY += 1
            elif token[i] == 'N':
                token[i] = 'N N'
                countN += 1
            else:
                warning("Non reconginzed alleles found!", token[i])
        temsnp.append(token)
    infile.close()
    return(temsnp)

# you can read a global without declaring it global, but to write a global, you need to declare it.
global countN
global countATCG
global countY
# snpmatrix = read_snp("AmesUSInbreds_AllZeaGBSv1.0_imputed_20130508_chr2.hmp.txt")    



def get_maf_missing(onesnp):
    snptable = Counter(onesnp)
    dict = {}
    if 'N N' in snptable:
        dict['missing'] = float(snptable['N N'])/4476
        del(snptable['N N'])
    else:
        dict['missing'] = 0
     
    if len(snptable) <= 1:
        dict['maf'] = 0
    elif len(snptable) == 2:
        ssnp = sorted(snptable.values())
        dict['maf'] = float(ssnp[0])/(ssnp[0] + ssnp[1])
        skey = sorted(snptable, key=snptable.__getitem__)
        dict['minor'] = skey[0].split()[0]
    elif len(snptable) == 3:
        mysnp = []
        for key in snptable:
            temsnp = key.split()
            if temsnp[0] != temsnp[1]
                heterkey = key
                hetervalue = snptable[key]
            mysnp.append(temsnp)
        if len(set(mysnp)) = 2:
            del(snptable[heterkey])
            ssnp = sorted(snptable.values())
            dict['maf'] = float(2*ssnp[0] + hetervalue)/(2*ssnp[0] + 2*hetervalue + 2*ssnp[1])
            skey = sorted(snptable, key=snptable.__getitem__)
            dict['minor'] = skey[0].split()[0]
            
        else:
            print(snptable)
          
    else:
        print(snptable)       
    return({})

def snp_merge(amesidx=, snpmatrix=):
    
    for i in range(len(amesidx)):
    
    
        if(amesidx[i][2] == 1):
            resout.append(snpmatrix[amesidx[i]-1])
        elif (amesidx[i][2] > 1):
            myidx = amesidx[i][1].split()
            for idx in myidx:
                mysnps.append(snpmstrix[idx -1])
            

def snpcondense(mysnp):
    snpset = set(mysnp)
    if 'N' in snpset:
        snpset.remove('N')
    
    snpset = list(snpset)
    if len(snpset) == 0:
        snpcall = "N"
    elif len(snpset) == 1:
        snpcall = snpset[0]
    elif len(snpset) == 2:
        snpcount = Counter(mysnp)
        if snpcount[snpset[0]] > snpcount[snpset[1]]:
            snpcall = snpset[0]
        else:
            snpcall = snpset[1]
    else:
        warning("multiple alleles found!")    
    return(snpcall)
        
        
    snpcount = Counter(mysnp)
    
    
    
    for snptype in snpset:
    


def removeN_set(mysnp):
    snpset = set(mysnp)
    if 'N' in snpset:
        snpset.remove('N')
    return snpset



##########################################################################################    
    





def readsnp_LofD(file_name, origin="hmp1"):

    infile = open(file_name, 'r')
    temp = []
    firstline = infile.readline()
    first = firstline.split()
    if first != header:
        warning("following this header:", header)
    secline = infile.readline()
    sec = secline.split()
    if(len(str(sec[4])) != 1):
        waring("SNP should code with one latter, like:'A'!")
    infile.close()

    infile = open(file_name, 'r')
    next(infile) # skip the header line
    for line in infile:
        tokens = line.split()
        tokens.append(origin)
        key = tokens[2] + "_" + tokens[3]
        temp.append({key:tokens})
    infile.close()
    return temp

def read_batch():
    print('Reading data from batch mode...')
    global hapmap1
    global hapmap2
    global rnaseq

    # Get the file names
    infile = open(args['infiles'], 'r')
    for line in infile:
        aline = line.split()
        if aline[0] == 'hapmap1':
            print("---- processing", os.path.basename(aline[1]))
            hapmap1 = hapmap1 + readsnp_LofD(aline[1], origin="hmp1")
        elif aline[0] == 'hapmap2':
            print("---- processing", os.path.basename(aline[1]))
            hapmap2 = hapmap2 + readsnp_LofD(aline[1], origin="hmp2")
        elif aline[0] == 'rnaseq':
            print("---- processing", os.path.basename(aline[1]))
            rnaseq = rnaseq + readsnp_LofD(aline[1], origin="rnaseq")
        else:
            warning("must indicate file type: hapmap1, hapmap2, rnaseq!")
    infile.close()

    input = hapmap1 + hapmap2 + rnaseq
    print(len(input), "SNPs were loaded!")
    return input

def read_chr():
    print('Reading HamMap1 data...')
    print('---- processing', os.path.basename(args['hmp1']))
    hapmap1 = readsnp_LofD(args['hmp1'], origin='hmp1')

    print('Reading HamMap2 data...')
    print('---- processing', os.path.basename(args['hmp2']))
    hapmap2 = readsnp_LofD(args['hmp2'], origin='hmp2')

    print('Reading RNA-seq data...')
    print('---- processing', os.path.basename(args['rnaseq']))
    rnaseq = readsnp_LofD(args['rnaseq'], origin='rnaseq')

    input = hapmap1 + hapmap2 + rnaseq
    print(len(input), "SNPs were loaded!")
    return input



##########################################################################################
"""
input is a list of dic
http://stackoverflow.com/questions/11358242/python-dictionary-with-same-keys
"""

def merge_snp():

    snpdic = {}
    for snp in input:

        ((x,y),) = snp.items()
        snpdic.setdefault(x,[]).append(y)

    for key, value in sorted(snpdic.items()):
        if len(value) == 1:
            unq[value[0][-1]].append(key)
            res.append(value[0])
        if len(value) == 2:
            tmp = comp2sets(akey=key, asnp=value)
            if tmp:
                res.append(tmp)
        if len(value) == 3:
            tmp = comp3sets(akey=key, asnp=value)
            if tmp:
                res.append(tmp)


def do_2snp_merge(htag='hmp12', asnp=['0']):
    snp1 = asnp[0]
    snp2 = asnp[1]
    merged_snp = asnp[0][0:4]
    for idx in range(4,31):
        temsnp = removeN_set([snp1[idx], snp2[idx]])
        temsnp = list(temsnp)
        if len(temsnp) == 2:
            merged_snp.append('N')
            aunmat[htag] += 1
        elif len(temsnp) == 0:
            merged_snp.append('N')
            amat[htag] += 1
        elif len(temsnp) == 1:
            merged_snp.append(temsnp[0])
            amat[htag] += 1
        else:
            warning(akey, "has multiple alleles:", temsnp)
    merged_snp.append(asnp[0][31])
    merged_snp.append(htag)
    return merged_snp

def comp2sets(akey='1_1000', asnp=['0']):

    snp1set = removeN_set(asnp[0][4:31])
    snp2set = removeN_set(asnp[1][4:31])

    sets = set([ asnp[0][-1], asnp[1][-1] ])

    merged_snp = []

    if snp1set == snp2set:
        if sets == set(["hmp1", "hmp2"]):
            lmat['hmp12'].append(akey)
            merged_snp = do_2snp_merge(htag='hmp12', asnp=asnp)
        elif sets == set(['hmp1', 'rnaseq']):
            lmat['hmp1rna'].append(akey)
            merged_snp = do_2snp_merge(htag='hmp1rna', asnp=asnp)
        elif sets == set(['hmp2', 'rnaseq']):
            lmat['hmp2rna'].append(akey)
            merged_snp = do_2snp_merge(htag='hmp2rna', asnp=asnp)
        else:
            warning(akey, "has set:", sets)

    else:
        if sets == set(["hmp1", "hmp2"]):
            lunmat['hmp12'].append(akey)
            merged_snp = do_2snp_merge(htag='hmp12', asnp=asnp)
        elif sets == set(['hmp1', 'rnaseq']):
            lunmat['hmp1rna'].append(akey)
            merged_snp = do_2snp_merge(htag='hmp1rna', asnp=asnp)
        elif sets == set(['hmp2', 'rnaseq']):
            lunmat['hmp2rna'].append(akey)
            merged_snp = do_2snp_merge(htag='hmp2rna', asnp=asnp)
        else:
            warning(akey, "has set:", sets)
        merged_snp = []

    return merged_snp



def removeN_set(mysnp):
    snpset = set(mysnp)
    if 'N' in snpset:
        snpset.remove('N')
    return snpset


def comp3sets(akey='1_1000', asnp=['0','1']):

    snp1set = removeN_set(asnp[0][4:31])
    snp2set = removeN_set(asnp[1][4:31])
    snp3set = removeN_set(asnp[2][4:31])

    mysnp = {}
    for tag in asnp:
        mysnp[tag[-1]] = tag

    merged_snp = []
    if snp1set == snp2set == snp3set:
        merged_snp = mysnp['hmp1'][0:4]
        lmat['hmp12rna'].append(akey)

        for idx in range(4,31):
            temsnp = removeN_set([mysnp['hmp1'][idx], mysnp['hmp2'][idx], mysnp['rnaseq'][idx]])
            temsnp = list(temsnp)
            if len(temsnp) == 2:
                aunmat['hmp12rna'] += 1
                merged_snp.append('N')
            elif len(temsnp) == 0:
                merged_snp.append('N')
            elif len(temsnp) == 1:
                amat['hmp12rna'] += 1
                merged_snp.append(temsnp[0])
            else:
                warning(akey, "has multiple alleles:", temsnp)

        merged_snp.append(mysnp['hmp1'][31])
        merged_snp.append("hmp12rna")
    else:
         lunmat['hmp12rna'].append(akey)

    return merged_snp

##########################################################################################

def writeMergedData():
	"""
	Write the merged data to a single file for use in downstream applications.

	Input: A string containing the path to write the file to.
	"""

	outfile = open(args['output'], 'w')

	outfile.write('\t'.join(header) + '\t' + 'source' + '\n')
	for snp in res:
		outfile.write('\t'.join(snp) + '\n')

	outfile.close()

def writeLog():
	"""
	Write a log file containing the initial SNP counts, the number of merged SNPs, and the final SNP count.

	Input: A string containing the path to write the log to.
	"""

	outfile = open(args['output'] + '.log', 'w')

	outfile.write('---Initial SNPs:---\n')
	outfile.write('Hapmap1: ' + str(len(hapmap1)) + '\nHapmap2: ' + str(len(hapmap2)) \
	 + '\nRNA-seq: ' + str(len(rnaseq)) + '\n')

	total = len(hapmap1) + len(hapmap2) + len(rnaseq)
	outfile.write('Total: ' + str(total) + '\n')

	merged_snps = total - len(res)
	outfile.write('---Merged SNPs: ' + str(merged_snps) + '\n')
	outfile.write('---Final SNPs: ' + str(len(res)) + '\n')
	outfile.write('\n')

	outfile.write('---Dataset unique loci---\n')
	outfile.write('-HapMap1 only: ' + str(len(unq['hmp1'])) + '\n')
	outfile.write('-HapMap2 only: ' + str(len(unq['hmp2'])) + '\n')
	outfile.write('-RNA-seq only: ' + str(len(unq['rnaseq'])) + '\n')

	outfile.write('---Dataset matched loci and alleles ---\n')
	outfile.write('-HapMap1 and HapMap2: ' + str(len(lmat['hmp12'])) + '\t' + str(amat['hmp12']) + '\n')
	outfile.write('-HapMap1 and RNA-seq: ' + str(len(lmat['hmp1rna'])) + '\t' + str(amat['hmp1rna']) + '\n')
	outfile.write('-HapMap2 and RNA-seq: ' + str(len(lmat['hmp2rna'])) + '\t' + str(amat['hmp2rna']) + '\n')
	outfile.write('-HapMap1 and HapMap2 and RNA-seq: ' + str(len(lmat['hmp12rna'])) + '\t' + str(amat['hmp12rna']) + '\n')

	outfile.write('---Dataset unmatched loci and alleles ---\n')
	outfile.write('-HapMap1 and HapMap2: ' + str(len(lunmat['hmp12'])) + '\t' + str(aunmat['hmp12']) + '\n')
	outfile.write('-HapMap1 and RNA-seq: ' + str(len(lunmat['hmp1rna'])) + '\t' + str(aunmat['hmp1rna']) + '\n')
	outfile.write('-HapMap2 and RNA-seq: ' + str(len(lunmat['hmp2rna'])) + '\t' + str(aunmat['hmp2rna']) + '\n')
	outfile.write('-HapMap1 and HapMap2 and RNA-seq: ' + str(len(lunmat['hmp12rna'])) + '\t' + str(aunmat['hmp12rna']) + '\n')

	outfile.close()


'''
main
'''
parser = get_parser()
args = vars(parser.parse_args())

if args['infiles'] or args['hmp1'] is None:
    print(version())

if args['path'] is not None:
    os.chdir(args['path'])

header = ['rs', 'alleles', 'chr', 'pos', 'B73', 'Z001', \
'Z002', 'Z003', 'Z004', 'Z005', 'Z006', 'Z007', 'Z008', 'Z009', \
'Z010', 'Z011', 'Z012', 'Z013', 'Z014', 'Z015', 'Z016', 'Z017', \
'Z018', 'Z019', 'Z020', 'Z021', 'Z022', 'Z023', 'Z024', 'Z025', \
'Z026', 'id' ]
hapmap1 = []
hapmap2 = []
rnaseq = []

#### read in the SNP data ####
if args['infiles'] is not None:
    input = read_batch()
else:
    input = read_chr()

#### output statistics ####
res = []
unq = {'hmp1':[], 'hmp2':[], 'rnaseq':[]}
lmat = {'hmp12':[], 'hmp1rna':[], 'hmp2rna':[], 'hmp12rna':[]}
lunmat = {'hmp12':[], 'hmp1rna':[], 'hmp2rna':[], 'hmp12rna':[]}
amat = {'hmp12':0, 'hmp1rna':0, 'hmp2rna':0, 'hmp12rna':0}
aunmat = {'hmp12':0, 'hmp1rna':0, 'hmp2rna':0, 'hmp12rna':0}


print('Merging the SNP data...')
merge_snp()

print('Writing merged data to file...')
writeMergedData()

print('Writing log file...')
writeLog()
