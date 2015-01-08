#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os

###############################
def readfile(infile_name="Land_hapguess_swith.out"):

    with open(infile_name, 'r') as infile:
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

test = readfile(infile_name="largedata/fphase/Parv_hapguess_switch.out")  

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
  
writePed(data=test, outfile="largedata/fphase/test.out")

### recycle this function from snpfrq
def get_loci_info(tokens):

    snptokens = tokens
    set0 = set(snptokens)
    if 'N' in set0:
        set0.remove('N')
    snpset = list(set0)
        
    info = {}
    if len(snpset) != 2:
      print "WARNING: ", "!=2 alleles !!!"
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


dtest = {'ind1':['A T C', 'C G A'], 
"ind2":['C T A', 'C G C'],
"ind3":['C T A', 'C G C']}

def transposeDofL(data):
  allhap = []
  for values in data.values():
    allhap.append(values[0].split())
    allhap.append(values[1].split())
  allsnp = map(list, zip(*allhap))
  
  snpinfo = []
  for i in range(len(allsnp)):
    snpinfo.append(get_loci_info(allsnp[i]))
  return snpinfo

res = transposeDofL(dtest)


#Genotype files for one chromosome in phased-EIGENSTRAT format (1 line per SNP): 
#Specified using parameters REFPOP1GENOFILE, REFPOP2GENOFILE
#Examples are CEUgenofile.22 and YRIgenofile.22
#Each line contains 2 columns per individual (00 or 01 or 10 or 11) 
#corresponding to the 1st and 2nd phased haplotype respectively.
#  e.g. for 4 samples and 2 SNPs we could have:
# 01110000
# 01000100

### determining minor and major alleles
def writeHAPMIX(data, snpinfo):
  
  for values in data.values():
    hap1 = values[0].split()
    hap2 = values[1].split()
  onesnp.extend([ hap1[i], hap2[i] ])
  
  for i in range(snpnum):
    onesnp = []
    
    allsnp.append(onesnp)
  return allsnp

test2 = writeHAPMIX(data=test, outfile="largedata/fphase/test.hapmix")




def readData(bedfile="largedata/IBD/allsnps_11m_IBD.bed", 
             gerpfile="largedata/SNP/allsnps_11m_gerpv2_tidy.csv", 
             dsffile="largedata/SNP/allsnps_11m.dsf5"):
  ### read.table
  snpibd = pd.read_table(bedfile, sep="\t", header= None)
  gerp = pd.read_csv(gerpfile)
  dsf7 = pd.read_table(dsffile, sep="\t")

  ### replace = gsub
  snpibd[3].replace("chr", "", regex=True, inplace=True)

  ### paste
  snpibd['snpid'] = snpibd[3].map(str) + "_" + snpibd[2].map(str)
  snpibd['ibdid'] = snpibd[3].map(str) + "_" + snpibd[4].map(str)

  ### merge SNPs with positive GERP score
  snpibd = pd.merge(snpibd, gerp[["snpid", "RS"]], on='snpid', sort=False, how='inner')
  snp0 = snpibd[snpibd['RS']>0]

  ###
  ibd0 = snp0.groupby('ibdid').size()
  #len(ibd0)
  #Out[16]: 85388

  ### GERP>0 merged with dsf7, dsf7 without B73 missing!
  dsf7 = dsf7[dsf7['B73'] != 'N']
  ibddsf = pd.merge(snp0, dsf7.iloc[:,np.r_[0, 7:19]], on="snpid", sort=False, how="inner")
  return ibddsf

### create a pedigree table
def getPed(ibddsf):
  names = ibddsf.columns.values[9:21]
  names = np.sort(names)
  ped = pd.DataFrame({'F1': np.random.randn(144),
                      'P1': np.repeat(names, 12, axis=0),
                      'P2': names.tolist() * 12})
  ped = ped[ped['P1'] < ped['P2']]
  ped.loc[:, 'F1'] = ped['P1'] + "x" + ped['P2']
  return ped
  
### compute GERP score in a given IBD block  
def ComputeOneGroup(onegroup):
  # gerp1a/gerp1d: gerp complementation
  # gerp2a/gerp2d: sum of gerp value
  gerp1a = gerp2a = gerp1d = gerp2d = 0
  p1 = "P1"
  p2 = "P2"
  for name, onesnp in onegroup.iterrows():
    # checking B73 reference
    if onesnp["B73"] != "N":
      b73 = onesnp["B73"]
        
      if onesnp[p1] == b73 and onesnp[p2] == b73:
        gerp1a = gerp1a + 2
        gerp2a = gerp2a + onesnp["RS"]*2
        gerp1d = gerp1d + 1
        gerp2d = gerp2d + onesnp["RS"]
      elif (onesnp[p1] != b73 and onesnp[p1] != "N") and onesnp[p2] == b73:
        gerp1a = gerp1a + 1
        gerp2a = gerp2a + onesnp["RS"]
        gerp1d = gerp1d + 1
        gerp2d = gerp2d + onesnp["RS"]
      elif onesnp[p1] == "N" and onesnp[p2] == b73:
        gerp1a = gerp1a + 1.5
        gerp2a = gerp2a + onesnp["RS"]*1.5
        gerp1d = gerp1d + 1
        gerp2d = gerp2d + onesnp["RS"]
          
      elif onesnp[p1] == b73 and (onesnp[p2] != b73 and onesnp[p2] != "N"):
        gerp1a = gerp1a + 1
        gerp2a = gerp2a + onesnp["RS"]
        gerp1d = gerp1d + 1
        gerp2d = gerp2d + onesnp["RS"]
      elif (onesnp[p1] != b73 and onesnp[p1] != "N") and (onesnp[p2] != b73 and onesnp[p2] != "N"):
        gerp1a = gerp1a + 0
        gerp2a = gerp2a + 0
        gerp1d = gerp1d + 0
        gerp2d = gerp2d + 0
      elif onesnp[p1] == "N" and (onesnp[p2] != b73 and onesnp[p2] != "N"):
        gerp1a = gerp1a + 0.5
        gerp2a = gerp2a + onesnp["RS"]*0.5
        gerp1d = gerp1d + 0.5
        gerp2d = gerp2d + onesnp["RS"]*0.5
       
      elif onesnp[p1] == b73 and onesnp[p2] == "N":
        gerp1a = gerp1a + 1.5
        gerp2a = gerp2a + onesnp["RS"]*1.5
        gerp1d = gerp1d + 1
        gerp2d = gerp2d + onesnp["RS"]
      elif (onesnp[p1] != b73 and onesnp[p1] != "N") and onesnp[p2] == "N":
        gerp1a = gerp1a + 0.5
        gerp2a = gerp2a + onesnp["RS"]*0.5
        gerp1d = gerp1d + 0.5
        gerp2d = gerp2d + onesnp["RS"]*0.5
      elif onesnp[p1] == "N" and onesnp[p2] == "N":
        gerp1a = gerp1a + 1
        gerp2a = gerp2a + onesnp["RS"]
        gerp1d = gerp1d + 1
        gerp2d = gerp2d + onesnp["RS"]
      else:
        warnings(onesnp["snpid"], "for", p1, p2, "have problem for additive imputation!")
    # for cases that B73==N
    else:
      gerp1a = gerp1a + 0
      gerp2a = gerp2a + 0
      gerp1d = gerp1d + 0
      gerp2d = gerp2d + 0
    gres = {'gerp1a': gerp1a, 'gerp2a': gerp2a, "gerp1d": gerp1d, "gerp2d": gerp2d}  
  return pd.Series(gres, name='metrics')
  

### Iterating through F1
def GetIBDgerp(ped, ibddsf):
  
  ### setup empty dataframe to collect results
  resa1 = resa2 = resd1 = resd2 = pd.DataFrame()
  
  for index, row in ped.iterrows():
    mydf = ibddsf[ ["ibdid", "B73", "RS", row["P1"], row["P2"]] ]
    mydf.columns = ["ibdid", "B73", "RS", 'P1', 'P2']
    
    print(">>> computing F1: [ ", row["F1"], " ]!")
    myres = mydf.groupby(['ibdid']).apply(ComputeOneGroup)
    #### concatenate the results
    tema1 = pd.DataFrame(myres, columns= ["gerp1a"])
    tema1.columns = [row["F1"]]
    resa1 = pd.concat([resa1, tema1], axis=1)
    
    tema2 = pd.DataFrame(myres, columns= ["gerp2a"])
    tema2.columns = [row["F1"]]
    resa2 = pd.concat([resa2, tema2], axis=1)
    
    temd1 = pd.DataFrame(myres, columns= ["gerp1d"])
    temd1.columns = [row["F1"]]
    resd1 = pd.concat([resd1, temd1], axis=1)
    
    temd2 = pd.DataFrame(myres, columns= ["gerp2d"])
    temd2.columns = [row["F1"]]
    resd2 = pd.concat([resd2, temd2], axis=1)
  
  hashres = {"gerpa1": resa1, "gerpa2":resa2, "gerpd1":resd1, "gerpd2":resd2}
  return hashres  

### write results
def writeRes(hashres, outbase="largedata/SNP/test"):
  #Apply operates on each row or column with the lambda function
  #axis = 0 -> act on columns, axis = 1 act on rows
  #x is a variable for the whole row or column
  #This line will scale minimum = 0 and maximum = 1 for each column
  #newrange = [-10, 10]
  #mfac = (newrange[1] - newrange[0])
  #change to (-10, 10)
  gerpa1 = hashres["gerpa1"]
  gerpa1 = gerpa1.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpa1 = np.round(gerpa1, 0)
  gerpa1 = gerpa1.transpose() 
  gerpa1.insert(0, "ibdid", gerpa1.index)
  #nm1 = list(gerpa1.columns.values)
  gerpa1.to_csv("_".join([outbase, "a1.gs"]), sep="\t", header=True, index=False, index_label=False)

  gerpa2 = hashres["gerpa2"]
  gerpa2 = gerpa2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 0)
  gerpa2 = np.round(gerpa2, 0)
  gerpa2 = gerpa2.transpose() 
  gerpa2.insert(0, "ibdid", gerpa2.index)
  gerpa2.to_csv("_".join([outbase, "a2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerpd1 = hashres["gerpd1"]
  gerpd1 = gerpd1.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpd1 = np.round(gerpd1, 0)
  gerpd1 = gerpd1.transpose()
  gerpd1.insert(0, "ibdid", gerpd1.index)
  gerpd1.to_csv("_".join([outbase, "d1.gs"]), sep="\t", header=True, index=False, index_label=False)

  gerpd2 = hashres["gerpd2"]
  gerpd2 = gerpd2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 0)
  gerpd2 = np.round(gerpd2, 0)
  gerpd2 = gerpd2.transpose() 
  gerpd2.insert(0, "ibdid", gerpd2.index)
  gerpd2.to_csv("_".join([outbase, "d2.gs"]), sep="\t", header=True, index=False, index_label=False)

def version():
    ver0 = """
    ##########################################################################################
    fastPHASE2other version 0.1
    Author: Jinliang Yang
    purpose: conversion fastPHASE to other formats
    --------------------------------

    updated: 1/7/2014
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

