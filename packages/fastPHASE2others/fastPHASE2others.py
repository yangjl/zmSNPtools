#!/usr/bin/env python

from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os

###### read the fastphase output to {'ind1':['A T C','C T G']}
def readfile(infile="largedata/fphase/fastphase_hapguess_swith.out"):

    with open(infile, 'r') as infile:
      data = {}
      nline = -9
      for line in infile:
        line = line.strip()
        if not line: continue
        if line.startswith("BEGIN GENOTYPES"):
          nline = 0
          continue
        if line.startswith("END GENOTYPES"):
          # define global variable `snpnum`
          snpnum = len(data[snpid][0].split())
          print ">>> loaded", len(data),  "individuals with", snpnum, "SNPs!"
          return data
        if nline >= 0 and nline % 3 == 0:
          data[line] = []
          snpid = line
          nline += 1
        elif nline >= 0 and nline % 3 != 0:
          data[snpid].append(line)
          nline += 1

# df = readfile(infile="land_94_chr10_hapguess_switch.out")

#A text file with no header line, and one line per individual with the following six fields:
#1.Individual's family ID ('FID')
#2.Individual's within-family ID ('IID'; cannot be '0')
#3.Within-family ID of father ('0' if father isn't in dataset)
#4.Within-family ID of mother ('0' if mother isn't in dataset)
#5.Sex code ('1' = male, '2' = female, '0' = unknown)
#6.Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

def writePed(data, outfile="largedata/fphase/test.out"):
  
  with open(outfile, "w") as outfile:
    for keys,values in data.items():
      outfile.write("\t".join([fid, keys, '0', '0', '0', '0'])) 
      hap1 = values[0].split()
      hap2 = values[1].split()
      for i in range(len(hap1)):
        outfile.write("\t" + "\t".join([hap1[i], hap2[i]]) )
      outfile.write("\n")
  
#writePed(data=test, outfile="largedata/fphase/test.out")

### recycled this function from snpfrq
### Note: missing data should change to N, or imputed
def get_loci_info(tokens):

    snptokens = tokens
    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    info = {}
    if len(snpset) > 2:
      print "WARNING: ", "!=2 alleles !!!", snpset
      sys.exit()
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

### transpose directionary of list, for example the test data
def transposeDofL(data):
  allhap = []
  for key in sorted(data.keys()):
    allhap.append(data[key][0].split())
    allhap.append(data[key][1].split())
  allsnp = map(list, zip(*allhap))
  
  snptable = {"allhap":allhap, "allsnp":allsnp}
  return snptable


#Genotype files for one chromosome in phased-EIGENSTRAT format (1 line per SNP): 
#Specified using parameters REFPOP1GENOFILE, REFPOP2GENOFILE
#Examples are CEUgenofile.22 and YRIgenofile.22
#Each line contains 2 columns per individual (00 or 01 or 10 or 11) 
#corresponding to the 1st and 2nd phased haplotype respectively.
#  e.g. for 4 samples and 2 SNPs we could have:
# 01110000
# 01000100
def getHAPMIX(snps, snpinfo):
  hapmixout = []
  if len(snps) != len(snpinfo):
    print "!!!error getHAPMIX, snps!=snpinfo", "\n"
    sys.exit()
  
  for i in range(len(snps)):
    hapmix = []
    minor = snpinfo[i]["minor"]
    
    for onesnp in snps[i]:
      if onesnp == minor:
        hapmix.append(1)
      else:
        hapmix.append(0)
    hapmixout.append(hapmix)  
  return hapmixout

#getHAPMIX(res["allsnp"], snpinfo)
### determining minor and major alleles
def writeHAPMIX(hapmixout,  outfile="test.out"):
  with open(outfile, "w") as outfile:
    for onesnp in hapmixout:
      outfile.write("".join(onesnp) + "\n")
  

#test2 = writeHAPMIX(data=test, outfile="largedata/fphase/test.hapmix")



def main():
  #### read data into dict
  data = readfile(infile="land_94_chr10_hapguess_switch.out")
  #### get two forms of the data
  #### res['allhap'] and res['allsnp']
  res = transposeDofL(data)
  #### get snp info
  snpinfo = []
  for onesnp in res["allsnp"]:
    snpinfo.append(get_loci_info(onesnp))
  #### prepare HAPMIX format
  hapmix = getHAPMIX(data, snpinfo)

  
  for i in range(len(snpinfo)):
    print snpinfo[i]['minor']
  
  for values in res.items():
    
  
  
  
  for i in range(snpnum):
    onesnp = []
    
    allsnp.append(onesnp)
  return allsnp




def version():
    ver0 = """
    ##########################################################################################
    fastPHASE2others version 0.2
    Author: Jinliang Yang
    purpose: convert fastPHASE to other formats
    --------------------------------

    updated: Jan. 21st, 2015
    ##########################################################################################
    """
    return ver0

#def warning(*objs):
#    print("WARNING: ", *objs, end='\n', file=sys.stderr)
#    sys.exit()

##########################################################################################
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
  parser.add_argument('-d','--ibd', help='IBD region in BED6 format', default="largedata/IBD/allsnps_11m_IBD.bed",  type=str)
  parser.add_argument('-s','--snp', help='founder SNP type in DSF5 format', default="largedata/SNP/allsnps_11m.dsf5",  type=str)
  parser.add_argument('-g','--gerp', help='gerp rates in csv format', default='largedata/SNP/allsnps_11m_gerpv2_tidy.csv', type=str)
  parser.add_argument('-o', '--output', help='base of the output file', default='gerpIBD_output', type=str)
  return parser
  #parser = get_parser()
  #parser.print_help()
    
def main():
  ### cal running time 
  st = timeit.default_timer()
  ### get the arguments
  parser = get_parser()
  args = vars(parser.parse_args())
  ### read data
  print(">>> reading data")
  if args['path'] is not None:
     os.chdir(args['path'])
  ibddsf = readData(bedfile=args['ibd'], 
                    gerpfile=args["gerp"], 
                    dsffile=args['snp'])
  print(">>> generating pedigree info")
  ### get ped info from names of idbdsf object
  ped = getPed(ibddsf)
  ### get IBM gerp looping through ped lines
  result = GetIBDgerp(ped, ibddsf)
  print(">>> writing results")
  writeRes(hashres=result, outbase=args['output'])
  
  ### get the end time
  et = timeit.default_timer()
  print(">>> [ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
  print("[[[ Job Done! ]]]")
 
if __name__ == "__main__":
  main()

