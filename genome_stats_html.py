#!/usr/bin/env python

## Defne Surujon
## 24 July 2016

## Summarize genome data - works best with Genbank files with more metadata

import os,sys
sys.path.append('/usr/lib/python2.7')
sys.path.append('/usr/lib/python2.7/dist-packages')
sys.path.append('/usr/lib/pymodules/python2.7')

import tempfile, shutil
from optparse import OptionParser
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

options = OptionParser(usage='%prog -i input -g genome_format -o outputdir -f htmlfile ',
                       description="Specify input genome file, its format (gbk or fasta) and output directory")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk or .fasta)")
options.add_option("-g","--genome_format",dest="gformat",
                   help="Input file format 'gbk' or 'fasta'")
options.add_option("-o","--outdir",dest="outputdir",
                   help="output directory")
options.add_option("-f","--htmlfile",dest="htmlfile",
                   help="output html")

#find string between two specified strings
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

#reverse complement of a string
def revcomp(s):
    sc=""
    for i in range(1,len(s)+1):
        if s[-i]=="T":sc+="A"
        if s[-i]=="A":sc+="T"
        if s[-i]=="G":sc+="C"
        if s[-i]=="C":sc+="G"
    return sc

#transpose a matrix
def transpose(m):
    mT=[[]for x in m[0]]
    for i in range(0,len(m[0])):
        for j in range(0,len(m)):
            mT[i].append(m[j][i])
    return mT
#This fasta reader assumes there is only one sequence in the genome file
def read_fasta(genome):
    f=open(genome)
    lines=f.readlines()
    f.close()
    genome_sequence=""
    genome_metadata={}
    for i in lines:
        if i[0]==">": genome_metadata["Header"]=i[1:].strip("\n")
        else: genome_sequence+=i.strip("\n")

    return genome_sequence, genome_metadata

#parse a genbank file. A sequence and metadata dictionary is returned. 
#the metadata dictionary has info about each one of the feature types
#specified. 
def read_gbk(genome):
    f=open(genome)
    lines=f.readlines()
    f.close()
    genome_sequence=""
    genome_metadata={}

    locusline=lines[0].split()
    genome_metadata["ACCN"]=locusline[1]
    if "circular" in locusline: genome_metadata["Type"]="circular"
    if "linear" in locusline: genome_metadata["Type"]="linear"

    isSequence=False
    isFeatures=False

    feature_types=["CDS","gene","rRNA","tRNA","misc_feature"]
    coords_dict={type:[[],[],[]] for type in feature_types}

    for line in lines:
        #grab source organism
        if "SOURCE" in line:
            genome_metadata["Source"]=" ".join(line.split()[1:])

        if "FEATURES" in line:
            isFeatures=True
            
        elif isFeatures:
            type_feature=line[:21].strip()
            
            if type_feature in feature_types:
                cmatrix,rc=get_feature_coords(line)
                coords_dict[type_feature][0]+=cmatrix[0]
                coords_dict[type_feature][1]+=cmatrix[1]
                coords_dict[type_feature][2]+=[revcomp for i in cmatrix[0]]

        #find the start of the sequence
        if "ORIGIN" in line:
            isSequence=True
            isFeatures=False
        #read the sequence until the ending "//"
        elif isSequence and line[0]!="/":

            chars=" 1234567890\n"
            for c in chars: line=line.replace(c,"")
            genome_sequence+=line
            #following is for python3 only
            #genome_sequence+=line.translate({ord(x):None for x in " 1234567890\n"})
    for feat in feature_types:
        genome_metadata[feat]=coords_dict[feat]

    return genome_sequence, genome_metadata
#if the location is formatted within a complement
#statement, retrieve the substring contained between
# "complement(" and ")"
def uncomplement(l):
    revcomp=False
    if l[:10]=="complement":
        revcomp=True
        lnew=l[11:-1]
    else:lnew=l
    return lnew,revcomp

#if the location is formatted as order(x...x)
#only get the first set of coordinates. 
def unorder(l):
    if l[:5]=="order":
        lnew=l[6:-1].split(",")[0]
    else:lnew=l
    return lnew

#if the location is formatted within a join statement
#return a list of coordinate pairs 
def unjoin(l,revcomp):
    if l[:4]=="join":
        lsubs=l[5:-1].split(",")
        lsubsnew=[[int(j) for j in uncomplement(i)[0].split("..")] for i in lsubs]
        if revcomp==False:revcomp=uncomplement(lsubs[0])[1]
    else:lsubsnew=[[int(i) for i in l.split("..")]]
    return lsubsnew,revcomp

#retrieve a list of coordinate pairs from the
#location field in the feature table
def get_feature_coords(l):
    revcomp=False
    #ambiguous cases are indicated by > or <
    #we simply ignore the ambiguity
    l=l.replace(">","")
    l=l.replace("<","")
    coords=l.split()[1]
    lc,revcomp=uncomplement(unorder(coords))
    lsubs,revcomp=unjoin(unorder(lc),revcomp)
    return transpose(lsubs),revcomp

#retrieve the nucleotide substrings from a set of coordinates
#and concatenate them into one sequence
def get_feature_nucs(seq,coords):
    feature_nucs=""
    for i in range(0,len(coords[0])):
        featstart=coords[0][i]-1
        featend=coords[1][i]
        featrc=coords[2][i]
        if featrc:feature_nucs+=revcomp(seq[featstart:featend])
        else:feature_nucs+=seq[featstart:featend]
    return feature_nucs

#Standard summary of nucleotide frequencies in a sequence. 
#provide a seqname for writing to file. 
def nuc_counts(seq, seqname):
	if len(seq)==0: return "<p>Length of "+seqname+" is 0, no nucleotide statistics are computed</p>"

	summary="<p>Nucleotide counts for "+seqname+":<p><table><tr><th>A</th><th>T</th><th>C</th><th>G</th><th>N</th></tr>"
	summary+="<tr><td>"+str(seq.count("A"))+"</td><td>"+str(seq.count("T"))+"</td><td>"+str(seq.count("C"))+"</td><td>"+str(seq.count("G"))+"</td><td>"+str(seq.count("N"))+"</td></tr></table>"

	pcN=str(100.*seq.count("N")/len(seq))
	pcGC=str(100.*(seq.count("C")+seq.count("G"))/len(seq))
	pcGCwoN=str(100.*(seq.count("C")+seq.count("G"))/(len(seq)-seq.count("N")))
	summary+="<p>%N's: "+pcN+"<br>%GC: "+pcGC+"<br>%GC (excluding N's): "+pcGCwoN+"<br></p>"

	return summary

#Make a histogram of [mydata], save it to [filename]
#specify the title of the plot to be written in the image
#file under [maintitle]
def histogram(mydata,filename,maintitle):
    fig,ax=plt.subplots()
    plt.hist(mydata,50)
    plt.xlabel("Length of sequence (nucleotide)")
    plt.ylabel("Count")
    plt.title(maintitle)
    #pl.show()
    plt.savefig(filename)


#generate a standard summary of a given feature type
#makes nucleotide frequency tables and histograms
def summarize_feature(sq,metadict,feature_name,outdir):
    html_segment="<h2>Feature: "+feature_name+"</h2>\n"
    feature_len=sum(metadict[feature_name][1])-sum(metadict[feature_name][0])
    html_segment+="<ul><li>Total length of "+feature_name+": "+str(feature_len)+"</li>"
    feature_nucs=get_feature_nucs(sq,metadict[feature_name])
    if len(metadict[feature_name][0])>0:
        feature_avg=feature_len/len(metadict[feature_name][0])
        html_segment+="<li>Average length of "+feature_name+": "+str(feature_avg)
    html_segment+="</li></ul>\n"+nuc_counts(feature_nucs,feature_name)
    if len(metadict[feature_name][0])>0:
        feature_lens=[metadict[feature_name][1][i]-metadict[feature_name][0][i] for i in range(0,len(metadict[feature_name][0]))]
        #Make a histogram for the coding sequences' lengths
        imgfile=outdir+"/"+feature_name+"hist.png"
        histogram(feature_lens,imgfile,"Feature: "+feature_name)
        html_segment+="<img src=\""+feature_name+"hist.png\">\n"
    return html_segment

def main():
    opts, args = options.parse_args()
    genomefile=opts.inputfile
    gformat=opts.gformat
    outputhtml=opts.htmlfile

    if not os.path.exists(opts.outputdir):
        os.makedirs(opts.outputdir)

    if gformat=="fasta":
        gnm, gnm_meta=read_fasta(genomefile)
    elif gformat=="gbk":
        gnm, gnm_meta=read_gbk(genomefile)
    else:
        print("Please specify an accepted genome file format('gbk' or 'fasta')")

    gnm=gnm.upper()
    gnm_meta["Acount"]=gnm.count("A")
    gnm_meta["Ccount"]=gnm.count("C")
    gnm_meta["Gcount"]=gnm.count("G")
    gnm_meta["Tcount"]=gnm.count("T")
    gnm_meta["Ncount"]=gnm.count("N")
    gnm_meta["Length"]=len(gnm)

    gnm_meta["%GC"]=100.*(gnm_meta["Gcount"]+gnm_meta["Ccount"])/gnm_meta["Length"]
    gnm_meta["%GCwoN"]=100.*(gnm_meta["Gcount"]+gnm_meta["Ccount"])/(gnm_meta["Length"]-gnm_meta["Ncount"])
    gnm_meta["%N"]=100.*(gnm_meta["Ncount"])/gnm_meta["Length"]

    #for i in gnm_meta: print(i,gnm_meta[i])

    f=open(outputhtml,"w")
    html_header="""
    <html><head><style>
    table, th, td {
        border: 1px solid black;
        border-collapse: collapse;
    }
    </style>
    </head>
    <body>
    """
    f.write(html_header)
    outputheader="<h1>Genome statistics for "+genomefile+"</h1>\n"
    f.write(outputheader+"\n")

    try:        #For FASTA
        f.write("<p>Genome header: "+gnm_meta["Header"]+"</p>\n")
    except KeyError:
        pass
    
    try:        #FOR GBK
        f.write("<p>Genome accession number: "+gnm_meta["ACCN"]+"<br>\n")
        f.write("Source organism: "+gnm_meta["Source"]+"<br>\n")
        f.write("Type: "+gnm_meta["Type"]+"</p>\n")
        
    except KeyError:
        pass
        
    f.write("<p>Length: "+str(gnm_meta["Length"])+" bp<br><br>\n")
    f.write("Nucleotide counts for the whole genome:</p>\n")
    f.write("<table><tr><th>A</th><th>T</th><th>C</th><th>G</th><th>N</th></tr>")
    f.write("<tr><td>"+str(gnm_meta["Acount"])+"</td><td>"+str(gnm_meta["Tcount"])+"</td><td>"+
			str(gnm_meta["Ccount"])+"</td><td>"+str(gnm_meta["Gcount"])+"</td><td>\n"+
			str(gnm_meta["Ncount"])+"</td></tr></table>\n")
    f.write("<p>%N's: "+str(gnm_meta["%N"])+"<br>")
    f.write("%GC: "+str(gnm_meta["%GC"])+"<br>")
    f.write("%GC (excluding N's): "+str(gnm_meta["%GCwoN"])+"</p>\n")

    try:        #for GBK
        feature_types=["CDS","gene","rRNA","tRNA","misc_feature"]
        #total number of features
        numfeats=0
        for i in feature_types:
            numfeats+=len(gnm_meta[i][0])

        f.write("<p>#Features: "+str(numfeats)+"</p>\n<ul>")
        for i in feature_types:
            f.write("<li>#"+i+": "+str(len(gnm_meta[i][0]))+"</li>")
        f.write("</ul>\n")

        for i in feature_types:
            f.write(summarize_feature(gnm,gnm_meta,i,opts.outputdir))

        
    except KeyError:
        pass
    
		
    f.write("</html>\n</body>")
    f.close()

if __name__ == '__main__':
    main()
