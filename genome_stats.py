## Defne Surujon
## 24 July 2016

## Parser for Genbank files


import os
from optparse import OptionParser


options = OptionParser(usage='%prog -i input -g genome_format -o output ',
                       description="Specify input genome file, its format (gbk or fasta) and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk or .fasta)")
options.add_option("-g","--genome_format",dest="gformat",
                   help="Input file format 'gbk' or 'fasta'")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.txt)")


#find string between two specified strings
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

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

    CDScoords=[[],[]]
    tRNAcoords=[[],[]]
    rRNAcoords=[[],[]]

    for line in lines:
        #grab source organism
        if "SOURCE" in line:
            genome_metadata["Source"]=" ".join(line.split()[1:])

        if "FEATURES" in line:
            isFeatures=True
            
        elif isFeatures:            #so far we are only inerested in 3 key feature types
            #CDS
            if "CDS" in line[:22]:
                cor1,cor2=get_feature_coords(line)
                CDScoords[0].append(cor1)
                CDScoords[1].append(cor2)
            #tRNA
            elif "tRNA" in line[:22]:
                cor1,cor2=get_feature_coords(line)
                tRNAcoords[0].append(cor1)
                tRNAcoords[1].append(cor2)
            #rRNA
            elif "rRNA" in line[:22]:
                cor1,cor2=get_feature_coords(line)
                tRNAcoords[0].append(cor1)
                tRNAcoords[1].append(cor2)

        #find the start of the sequence
        if "ORIGIN" in line:
            isSequence=True
            isFeatures=False
        #read the sequence until the ending "//"
        elif isSequence and line[0]!="/":
            genome_sequence+=line.translate({ord(x):None for x in " 1234567890\n"})

    genome_metadata["CDS"]=CDScoords
    genome_metadata["tRNA"]=tRNAcoords
    genome_metadata["rRNA"]=rRNAcoords
    
    return genome_sequence, genome_metadata

def get_feature_coords(l):
    coords=l.split()[1]
    #print (coords)
    if "complement" in coords:
        c2=find_between(coords,"..",")")
        c1=find_between(coords,"(","..")
    else:
        c1=coords.split("..")[0]
        c2=coords.split("..")[1]
    return int(c1), int(c2)
    
def main():
    opts, args = options.parse_args()
    genomefile=opts.inputfile
    gformat=opts.gformat
    outputtxt=opts.outputfile

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
    gnm_meta["%N"]=100.*(gnm_meta["Ncount"])/gnm_meta["Length"]

    #for i in gnm_meta: print(i,gnm_meta[i])

    f=open(outputtxt,"w")
    outputheader="Genome statistics for "+genomefile
    f.write(outputheader+"\n")
    f.write("="*len(outputheader)+"\n")
    
    try:        #For FASTA
        f.write("Genome header: "+gnm_meta["Header"]+"\n")
    except KeyError:
        pass
    
    try:        #FOR GBK
        f.write("Genome accession number: "+gnm_meta["ACCN"]+"\n")
        f.write("Source organism: "+gnm_meta["Source"]+"\n")
        f.write("Type: "+gnm_meta["Type"]+"\n")
        
    except KeyError:
        pass
        
    f.write("Length: "+str(gnm_meta["Length"])+" bp\n\n")
    f.write("Nucleotide counts for the whole genome:\n")
    f.write("#A\t#T\t#C\t#G\t#N\n")
    f.write(str(gnm_meta["Acount"])+"\t"+str(gnm_meta["Tcount"])+"\t"+
            str(gnm_meta["Ccount"])+"\t"+str(gnm_meta["Gcount"])+"\t"+
            str(gnm_meta["Ncount"])+"\n")
    f.write("%N's: "+str(gnm_meta["%N"])+"\n")
    f.write("%GC: :"+str(gnm_meta["%GC"])+"\n\n")

    try:        #for GBK
        #total number of features
        numfeats=len(gnm_meta["CDS"][0])+len(gnm_meta["tRNA"][0])+len(gnm_meta["rRNA"][0])

        f.write("#Features: "+str(numfeats)+"\n")
        f.write("\t#CDS: "+str(len(gnm_meta["CDS"][0]))+"\n")
        f.write("\t#tRNA: "+str(len(gnm_meta["tRNA"][0]))+"\n")
        f.write("\t#rRNA: "+str(len(gnm_meta["rRNA"][0]))+"\n")

        #lengths of features
        CDSlen=sum(gnm_meta["CDS"][1])-sum(gnm_meta["CDS"][0])
        f.write("Total length of CDS: "+str(CDSlen)+"\n")
        if len(gnm_meta["CDS"][0])>0:
            CDSavg=CDSlen/len(gnm_meta["CDS"][0])
            f.write("Average length of CDS: "+str(CDSavg)+"\n")

        tRNAlen=sum(gnm_meta["tRNA"][1])-sum(gnm_meta["tRNA"][0])
        f.write("Total length of tRNA: "+str(tRNAlen)+"\n")
        if len(gnm_meta["tRNA"][0])>0:
            tRNAavg=tRNAlen/len(gnm_meta["tRNA"][0])
            f.write("Average length of tRNA: "+str(tRNAavg)+"\n")

        rRNAlen=sum(gnm_meta["rRNA"][1])-sum(gnm_meta["rRNA"][0])
        f.write("Total length of rRNA: "+str(rRNAlen)+"\n")
        if len(gnm_meta["rRNA"][0])>0:
            rRNAavg=rRNAlen/len(gnm_meta["rRNA"][0])
            f.write("Average length of rRNA: "+str(rRNAavg)+"\n")
        
    except KeyError:
        pass


    f.close()

if __name__ == '__main__':
    main()
