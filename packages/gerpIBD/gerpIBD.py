#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import pandas as pd
import numpy as np
import sys
import argparse
import textwrap
import timeit
import os

def readData(bedfile="largedata/IBD/allsnps_11m_IBD.bed", 
             gerpfile="largedata/SNP/gerpv2_b0_real.csv", 
             dsffile="largedata/SNP/allsnps_11m.dsf5",
             hfile="largedata/snpeff/gy_h.txt",
             gerppositive="positive"):
  ### read.table
  snpibd = pd.read_table(bedfile, sep="\t", header= None)
  gerp = pd.read_csv(gerpfile)
  dsf7 = pd.read_table(dsffile, sep="\t")
  h = pd.read_table(hfile, sep="\t")
  
  ### merge gerp with h
  gerph = pd.merge(gerp, h, on='snpid', sort=False, how="inner")
  
  ### replace = gsub
  snpibd[3].replace("chr", "", regex=True, inplace=True)

  ### paste
  snpibd['snpid'] = snpibd[3].map(str) + "_" + snpibd[2].map(str)
  snpibd['ibdid'] = snpibd[3].map(str) + "_" + snpibd[4].map(str)

  ### merge SNPs with positive GERP H score
  snpibd = pd.merge(snpibd, gerph[["snpid", "RS", "h"]], on='snpid', sort=False, how='inner')
  if gerppositive == "positive":
    snp0 = snpibd[snpibd['RS']>0]
  else:
    snp0 = snpibd

  ### check how many SNP per IBD block ###
  ibd0 = snp0.groupby('ibdid').size()
  #len(ibd0)
  #Out[16]: 85388

  ### GERP>0 merged with dsf7, dsf7 without B73 missing!
  # dsf7 = dsf7[dsf7['B73'] != 'N']
  ibddsf = pd.merge(snp0, dsf7.iloc[:,np.r_[0, 3, 7:19]], on="snpid", sort=False, how="inner")
  return ibddsf

### create a pedigree table
def getPed(ibddsf):
  names = ibddsf.columns.values[11:23]
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
  gerp2a = gerp2d = gerp2h = 0
  gerpa2b = gerpab2 = 0
  p1 = "P1"
  p2 = "P2"
  for name, onesnp in onegroup.iterrows():
    # checking B73 reference
    if onesnp["major"] != "N":
      b73 = onesnp["major"]
        
      if onesnp[p1] == b73 and onesnp[p2] == b73:
        gerp2a = gerp2a + onesnp["RS"]*2
        gerp2d = gerp2d + onesnp["RS"]
        gerp2h = gerp2h + onesnp["RS"]
        gerpa2b = gerpa2b + onesnp["RS"]*2
        gerpab2 = gerpab2 + onesnp["RS"]*2
      elif (onesnp[p1] != b73 and onesnp[p1] != "N") and onesnp[p2] == b73: 
        gerp2a = gerp2a + onesnp["RS"]
        gerp2d = gerp2d + onesnp["RS"]
        gerp2h = gerp2h + onesnp["RS"]*(1+onesnp["h"])
        gerpa2b = gerpa2b + onesnp["RS"]*(1+onesnp["h"])*2/3
        gerpab2 = gerpab2 + onesnp["RS"]*(1+onesnp["h"])*4/3
      elif onesnp[p1] == "N" and onesnp[p2] == b73:
        gerp2a = gerp2a + onesnp["RS"]*1.5
        gerp2d = gerp2d + onesnp["RS"]
        gerp2h = gerp2h + onesnp["RS"]*(1.5+onesnp["h"]/2)
        gerpa2b = gerpa2b + onesnp["RS"]*(4/3+onesnp["h"]*1/3)
        gerpab2 = gerpab2 + onesnp["RS"]*(5/3+onesnp["h"]*2/3)
      elif onesnp[p1] == b73 and (onesnp[p2] != b73 and onesnp[p2] != "N"):
        gerp2a = gerp2a + onesnp["RS"]
        gerp2d = gerp2d + onesnp["RS"]
        gerp2h = gerp2h + onesnp["RS"]*(1+onesnp["h"])
        gerpa2b = gerpa2b + onesnp["RS"]*(1+onesnp["h"])*2/3
        gerpab2 = gerpab2 + onesnp["RS"]*(1+onesnp["h"])*4/3
      elif (onesnp[p1] != b73 and onesnp[p1] != "N") and (onesnp[p2] != b73 and onesnp[p2] != "N"):
        gerp2a = gerp2a + 0
        gerp2d = gerp2d + 0
        gerp2h = gerp2h + 0
        gerpa2b = gerpa2b + 0
        gerpab2 = gerpab2 + 0
      elif onesnp[p1] == "N" and (onesnp[p2] != b73 and onesnp[p2] != "N"):
        gerp2a = gerp2a + onesnp["RS"]*0.5
        gerp2d = gerp2d + onesnp["RS"]*0.5
        gerp2h = gerp2h + onesnp["RS"]*(1+onesnp["h"])/2
        gerpa2b = gerpa2b + onesnp["RS"]*(1+onesnp["h"])/3
        gerpab2 = gerpab2 + onesnp["RS"]*(1+onesnp["h"])/3*2
      elif onesnp[p1] == b73 and onesnp[p2] == "N":
        gerp2a = gerp2a + onesnp["RS"]*1.5
        gerp2d = gerp2d + onesnp["RS"]
        gerp2h = gerp2h + onesnp["RS"]*(1.5+onesnp["h"]/2)
        gerpa2b = gerpa2b + onesnp["RS"]*(4+onesnp["h"])/3
        gerpab2 = gerpab2 + onesnp["RS"]*(5+onesnp["h"]*2)/3
      elif (onesnp[p1] != b73 and onesnp[p1] != "N") and onesnp[p2] == "N":
        gerp2a = gerp2a + onesnp["RS"]*0.5
        gerp2d = gerp2d + onesnp["RS"]*0.5
        gerp2h = gerp2h + onesnp["RS"]*onesnp["h"]*0.5
        gerpa2b = gerpa2b + onesnp["RS"]*onesnp["h"]*0.5*1/3
        gerpab2 = gerpab2 + onesnp["RS"]*onesnp["h"]*0.5*2/3
      elif onesnp[p1] == "N" and onesnp[p2] == "N":
        gerp2a = gerp2a + onesnp["RS"]
        gerp2d = gerp2d + onesnp["RS"]*2/3
        gerp2h = gerp2h + onesnp["RS"]*(1+onesnp["h"]/3)
        gerpa2b = gerpa2b + onesnp["RS"]*(8/9+2/9*onesnp["h"])
        gerpab2 = gerpab2 + onesnp["RS"]*(10/9+4/9*onesnp["h"])
      else:
        warnings(onesnp["snpid"], "for", p1, p2, "have problem for additive imputation!")
    # for cases that B73==N
    else:
      gerp2a = gerp2a + 0
      gerp2d = gerp2d + 0
      gerp2h = gerp2h + 0
      gerpa2b = gerpa2b + 0
      gerpab2 = gerpab2 + 0
    gres = {'gerp2a': gerp2a, "gerp2d": gerp2d, "gerp2h": gerp2h, "gerpa2b":gerpa2b, "gerpab2":gerpab2}  
  return pd.Series(gres, name='metrics')
  

### Iterating through F1
def GetIBDgerp(ped, ibddsf):
  
  ### setup empty dataframe to collect results
  resa2 = resd2 = resh2 = a2b = ab2 = pd.DataFrame()
  
  for index, row in ped.iterrows():
    mydf = ibddsf[ ["ibdid", "major", "RS", "h", row["P1"], row["P2"]] ]
    mydf.columns = ["ibdid", "major", "RS", "h", 'P1', 'P2']
    
    print(">>> computing F1: [ ", row["F1"], " ]!")
    myres = mydf.groupby(['ibdid']).apply(ComputeOneGroup)
    #### concatenate the results
    
    tema2 = pd.DataFrame(myres, columns= ["gerp2a"])
    tema2.columns = [row["F1"]]
    resa2 = pd.concat([resa2, tema2], axis=1)
    
    temd2 = pd.DataFrame(myres, columns= ["gerp2d"])
    temd2.columns = [row["F1"]]
    resd2 = pd.concat([resd2, temd2], axis=1)
    
    temh2 = pd.DataFrame(myres, columns= ["gerp2h"])
    temh2.columns = [row["F1"]]
    resh2 = pd.concat([resh2, temh2], axis=1)
    
    tema2b = pd.DataFrame(myres, columns= ["gerpa2b"])
    tema2b.columns = [row["F1"]]
    a2b = pd.concat([a2b, tema2b], axis=1)
    
    temab2 = pd.DataFrame(myres, columns= ["gerpab2"])
    temab2.columns = [row["F1"]]
    ab2 = pd.concat([ab2, temab2], axis=1)
  
  hashres = {"gerpa2":resa2, "gerpd2":resd2, "gerph2":resh2, "a2b":a2b, "ab2":ab2}
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
  
  gerpa2 = hashres["gerpa2"]
  #gerpa2 = gerpa2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpa2 = np.round(gerpa2, 0)
  gerpa2 = gerpa2.transpose() 
  gerpa2.insert(0, "ibdid", gerpa2.index)
  gerpa2.to_csv("_".join([outbase, "a2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerpd2 = hashres["gerpd2"]
  #gerpd2 = gerpd2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpd2 = np.round(gerpd2, 0)
  gerpd2 = gerpd2.transpose() 
  gerpd2.insert(0, "ibdid", gerpd2.index)
  gerpd2.to_csv("_".join([outbase, "d2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerph2 = hashres["gerph2"]
  #gerph2 = gerph2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerph2 = np.round(gerph2, 0)
  gerph2 = gerph2.transpose() 
  gerph2.insert(0, "ibdid", gerph2.index)
  gerph2.to_csv("_".join([outbase, "h2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerpa2b = hashres["a2b"]
  #gerpa2b = gerpa2b.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpa2b = np.round(gerpa2b, 0)
  gerpa2b = gerpa2b.transpose() 
  gerpa2b.insert(0, "ibdid", gerpa2b.index)
  gerpa2b.to_csv("_".join([outbase, "a2b.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerpab2 = hashres["ab2"]
  #gerpab2 = gerpab2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpab2 = np.round(gerpab2, 0)
  gerpab2 = gerpab2.transpose() 
  gerpab2.insert(0, "ibdid", gerpab2.index)
  gerpab2.to_csv("_".join([outbase, "ab2.gs"]), sep="\t", header=True, index=False, index_label=False)

### write results
def write_adk_only(hashres, outbase="largedata/SNP/test"):
  #Apply operates on each row or column with the lambda function
  #axis = 0 -> act on columns, axis = 1 act on rows
  #x is a variable for the whole row or column
  #This line will scale minimum = 0 and maximum = 1 for each column
  #newrange = [-10, 10]
  #mfac = (newrange[1] - newrange[0])
  #change to (-10, 10)
  
  gerpa2 = hashres["gerpa2"]
  #gerpa2 = gerpa2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpa2 = np.round(gerpa2, 0)
  gerpa2 = gerpa2.transpose() 
  gerpa2.insert(0, "ibdid", gerpa2.index)
  gerpa2.to_csv("_".join([outbase, "a2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerpd2 = hashres["gerpd2"]
  #gerpd2 = gerpd2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerpd2 = np.round(gerpd2, 0)
  gerpd2 = gerpd2.transpose() 
  gerpd2.insert(0, "ibdid", gerpd2.index)
  gerpd2.to_csv("_".join([outbase, "d2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
  gerph2 = hashres["gerph2"]
  #gerph2 = gerph2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerph2 = np.round(gerph2, 0)
  gerph2 = gerph2.transpose() 
  gerph2.insert(0, "ibdid", gerph2.index)
  gerph2.to_csv("_".join([outbase, "h2.gs"]), sep="\t", header=True, index=False, index_label=False)

### write results
def write_k_only(hashres, outbase="largedata/SNP/test"):
  #Apply operates on each row or column with the lambda function
  #axis = 0 -> act on columns, axis = 1 act on rows
  #x is a variable for the whole row or column
  #This line will scale minimum = 0 and maximum = 1 for each column
  #newrange = [-10, 10]
  #mfac = (newrange[1] - newrange[0])
  #change to (-10, 10)
  
  gerph2 = hashres["gerph2"]
  #gerph2 = gerph2.apply(lambda x:-10+(x.astype(float) - min(x))/(max(x)-min(x))*20, axis = 1)
  gerph2 = np.round(gerph2, 0)
  gerph2 = gerph2.transpose() 
  gerph2.insert(0, "ibdid", gerph2.index)
  gerph2.to_csv("_".join([outbase, "h2.gs"]), sep="\t", header=True, index=False, index_label=False)
  
def version():
    ver0 = """
    ##########################################################################################
    gerpIBD version 1.1
    Author: Jinliang Yang
    purpose: compute the accumulative GERP rate in an IBD region
    --------------------------------
    
    updated: 04/06/2016, change output mode.
    updated: 04/06/2016, change major rather than B73 allele as beneficial allele!
    updated: 09/25/2015, add argument outtype; removed a1 and d1, no normalization!
    updated: 09/20/2015, imputation for triplotype, a2b and ab2
    updated: 09/08/2015, incomplete dominance
    updated: 2/22/2014, do negative gerp
    ##########################################################################################
    """
    return ver0

def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

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
  parser.add_argument('-s','--snp', help='founder SNP type in DSF5 format, Note: major=should be beneficial allele.', 
                      default="largedata/SNP/allsnps_11m.dsf5",  type=str)
  parser.add_argument('-g','--gerp', help='gerp rates in csv format', default='largedata/SNP/allsnps_11m_gerpv2_tidy.csv', type=str)
  parser.add_argument('-f','--dofd', help='degree of dominance', default='largedata/snpeff/gy_h.txt', type=str)
  parser.add_argument('-n', '--num', help='Only use positive numbers of GERP', default='positive', type=str)
  parser.add_argument('-o', '--output', help='base of the output file', default='gerpIBD_output', type=str)
  parser.add_argument('-t', '--outtype', help='all, output all; adk output add, dom and k; k output k only!', default= 'all', type=str)
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
  print("###>>> reading data ...")
  if args['path'] is not None:
      os.chdir(args['path'])
  ibddsf = readData(bedfile=args['ibd'], 
                    gerpfile=args["gerp"], 
                    dsffile=args['snp'],
                    hfile=args['dofd'],
                    gerppositive=args['num'])
  print("###>>> generating pedigree info ...")
  ### get ped info from names of idbdsf object
  ped = getPed(ibddsf)
  ### get IBM gerp looping through ped lines
  result = GetIBDgerp(ped, ibddsf)
  print("###>>> writing results ...")
  if(args['outtype'] == "all"):
      writeRes(hashres=result, outbase=args['output'])
  elif(args['outtype'] == "adk"):
      write_adk_only(hashres=result, outbase=args['output'])
  elif(args['outtype'] == "k"):
      write_k_only(hashres=result, outbase=args['output'])
      
  ### get the end time
  et = timeit.default_timer()
  print("###>>> [ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
  print("###>>> [[[ Job Done! ]]]")
 
if __name__ == "__main__":
  main()

