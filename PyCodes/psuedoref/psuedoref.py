## This script takes reference fasta and a SNP vcf file and creates strain-specific fasta by substituting for SNPs. 
## Use of strain-specific genome removes allelic bias during read alignment. 
## This script is written for creating strain-specific reference file using Sanger MGP vcf file consiting of 18 strains. You can give the index number of strain for which SNPs need to be substituted. 

import sys,os,re,fileinput

Argument = []
Argument = sys.argv[1:]

filename = Argument[0] #Reference fasta
snpfile = Argument[1] #SNP file
strain = Argument[2] 

def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in xrange(0, len(string), every))


def substitute_snp(header,genome,SNP):
    name = ""
name = header
genomearray = []
genomearray = list(genome.rstrip(" "))

for snp in SNP:
    #print snp
    #print genomearray[int(snp)-1]
    #print SNP[snp][0]
    
    if (genomearray[int(snp)-1]).lower() == (SNP[snp][0]).lower():
    genomearray[int(snp)-1] = SNP[snp][1]
return  header+"\n"+insert_newlines("".join(genomearray))+"\n"


def nosnp_format(header,genome):
    name = ""
name = header
genomearray = []
genomearray = list(genome.rstrip(" "))
return  header+"\n"+insert_newlines("".join(genomearray))+"\n"


SNP = {}

def trim(list):
    listn = []
listn = list
newlist = []
for a in listn:
    newlist.append((a.rstrip()).lstrip())
return newlist


for line in fileinput.input([snpfile]):
    if line.startswith("##") or not line.strip():
    continue
else:
    info1 = []
info1 = line.rstrip("\n").split("\t")

info = []
info = trim(info1)
if info[strain] != "." and info[strain].startswith("1/1:1"):
    if "," not in info[4]:
    if info[0] not in SNP:
    SNP[info[0]] = {}
SNP[info[0]][info[1]] = []
SNP[info[0]][info[1]].append(info[3])
SNP[info[0]][info[1]].append(info[4])
else:
    SNP[info[0]][info[1]] = []
SNP[info[0]][info[1]].append(info[3])
SNP[info[0]][info[1]].append(info[4])

#print SNP

mousefile = open("Mus_musculus_SNPs.fa","w") #output file 

genomefile = open(filename, 'r').read() #reading whole file in one string
genomeb = []
dbafasta = ""

genomeb = genomefile.split(">")

for chr in genomeb[1:]:
    fasta = ""
dbfasta = ""
header = ""
genome = []
genome = chr.split("\n")

header = ">"+str(genome[0])
print genome[0]

if (genome[0].rstrip(" ")).lstrip("chr") in SNP:
    fasta = "".join(genome[1:])
dbafasta = substitute_snp(header,fasta,SNP[(genome[0].rstrip(" ")).lstrip("chr")])
mousefile.write(str(dbafasta))
else:
    fasta = "".join(genome[1:])
dbafasta = nosnp_format(header,fasta)
mousefile.write(str(dbafasta))

mousefile.close()