#!/usr/bin/env python
# File created on 21 Jun 2016

__author__ = "Defne Surujon"
__credits__=["Defne Surujon","Alexander Farrell","Federico Rosconi","Jon Anthony","Tim van Opijnen"]
__version__ = "0.0.1"
__maintainer__="Defne Surujon"
__email__="surujon@bc.edu"


import os
import subprocess
import datetime
from optparse import OptionParser

script_info = {}
script_info['brief_description'] = """Assemble a collection of contigs into a complete genome"""
script_info['script_description'] = \
"""This script assembles a complete genome by putting together a number of contigs
provided in a single fasta file. The contigs are compared to 20 complete Streptococcus
pneumoniae genomes, and the genome with the highest number of nucleotide matches is 
selected as the reference. The contigs are then placed on this reference genome (either
in the forward orientation or in reverse complement, as necessary)and spaced accordingly 
with an appropriate number of "N"s. The output is a single fasta file of the complete
genome.
"""

options = OptionParser(usage='%prog -i inputfile',
                       description="Specify sequence file and generate a complete genbank genome")
options.add_option("-i","--infile",dest="inputfile",
                   help="Input file of contigs(.fasta)")
options.add_option("-l","--logfile",dest="logfile",
			help="Log file (.txt)")

#lengths of genomes in the reference database
genome_lens={"NC_003028":2160842,"NC_014498":2240045,"NC_018630":2064154,"NC_003098":2038615,
	"NC_012469":2112148,"NC_010380":2245615,"NC_012468":2184682,"NC_012466":2120234,
	"NC_012467":2111882,"NC_011072":2078953,"NC_010582":2209198,"NC_014251":2088772,
	"NC_011900":2221315,"NC_014494":2130580,"NC_017593":2093317,"NC_017591":2142122,
	"NC_018594":2129934,"NC_017769":2150813,"NC_008533":2046115,"NC_017592":2036867}

#calculate the number of Ns coming after contig i
def numberN(order,c,cs,ce,gs,ge,f,i):
    n=0
    if f[i]==1 and f[i+1]==1:
        n=(gs[i+1]-ge[i])-cs[i+1]-(len(c[order[i]])-ce[i])
    elif f[i]==0 and f[i+1]==1:
        n=(gs[i+1]-gs[i])-cs[i]-cs[i+1]
    elif f[i]==1 and f[i+1]==0:
        n=(ge[i+1]-ge[i])-(len(c[order[i]])-ce[i])-(len(c[order[i+1]])-ce[i+1])
    elif f[i]==0 and f[i+1]==0:
        n=(ge[i+1]-gs[i])-cs[i]-(len(c[order[i+1]])-ce[i+1])

    return n

#reverse complement of a string
def revcomp(s):
    sc=""
    for i in range(1,len(s)+1):
        if s[-i]=="T":sc+="A"
        if s[-i]=="A":sc+="T"
        if s[-i]=="G":sc+="C"
        if s[-i]=="C":sc+="G"
    return sc

#read fasta file, return list of contigs
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles=[]
    myseq=""
    for i in lines[1:]:
        if i[0]==">":
            titles.append(myseq)
            myseq=""
        else:
            myseq+=(i.upper()).replace("\n","")
    #don't forget to grab the last seq
    titles.append(myseq)
    return titles

#write a new fasta file with the appropriate combination of contigs
def write_contigs(c,order,cstart,cend,gstart,gend,rev,filename,gen):
    thisfile=filename+"_COMPLETE.fasta"
    f=open(thisfile,'w')

    f.write(">"+filename+".complete\n")
    #leading N's:
    if rev[0]==1:
        f.write("N"*(gstart[0]-cstart[0]))
    else:
        f.write("N"*(gstart[0]+cend[0]-len(c[order[0]])))
    nums=0
    for i in range(0,len(c)):
        
        if rev[i]==1:
            numslast=(gen-gend[i])-(len(c[order[i]])-cend[i])
            f.write(c[order[i]])
        elif rev[i]==0:
            numslast=(gen-gend[i])-cstart[i]
            f.write(revcomp(c[order[i]]))
        
        if i==len(c)-1:nums=numslast
        else:nums=numberN(order,c,cstart,cend,gstart,gend,rev,i)
        f.write("N"*max(nums,25))

    f.close()

#read the preferences file, return two lists
def read_pref(preffile):
    f2=open(preffile)
    f2lines=f2.readlines()
    f2.close()
    o=[int(x)-1 for x in f2lines[0].split()]
    cstart=[int(x) for x in f2lines[1].split()]
    cend=[int(x) for x in f2lines[2].split()]
    gstart=[int(x) for x in f2lines[3].split()]
    gend=[int(x) for x in f2lines[4].split()]
    r=[int(x) for x in f2lines[5].split()]

    return o,cstart,cend,gstart,gend,r

def process_fasta(inputfile, preferences,genome):
    allcontigs=read_sequences(inputfile)
    ordering, cstart,cend,gstart,gend,rc=read_pref(preferences)
    write_contigs(allcontigs,ordering,cstart,cend,gstart,gend,rc,inputfile,genome)
	
def main():
	opts, args = options.parse_args()
	filename=str(opts.inputfile)
	dboutput=filename+"_results.txt"
	
	log=open(str(opts.logfile),"w")
	
	log.write("Starting BLAST search at "+str(datetime.datetime.now())+"\n")

	blasttask2=["blastn", "-query", filename, "-db", "/data1/BioData/BlastDBs/StrepStrains/strep-strains",
	"-out", dboutput, "-strand", "both", "-word_size", "8", "-evalue", "1E-20", "-outfmt",
	"10 qseqid qstart qend evalue sseqid sstart send"]
	subprocess.call(blasttask2, shell=False)
	
	log.write("BLAST run ended at "+str(datetime.datetime.now())+"\n")
	log.write("Rearranging contigs\n")
	
	subprocess.call(["Rscript", "get_longest_match_loci.R", dboutput,str(opts.logfile)])

	ftemp=open(str(opts.logfile),"r")
	templines=ftemp.readlines()
	ftemp.close()
	refgenome=templines[-3].split("|")[-2]

	process_fasta(filename,dboutput+"_positions.txt",genome_lens[refgenome[:-2]])

	log.write("Run ended at "+str(datetime.datetime.now())+"\n")
	log.close()

if __name__ == '__main__':
    main()