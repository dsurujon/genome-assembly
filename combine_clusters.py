import os
from optparse import OptionParser
import csv
import sys

options = OptionParser(usage='%prog -i input_dir -o output ', description = "Specify input cluster directory (-i) and output CSV file (-o)")

options.add_option("-i","--indir",dest="inputdir",help="Input directory")
options.add_option("-o","--outfile",dest="outputfile",help="Output CSV file (.csv)")
options.add_option("-d","--debug", dest="isDebug", action="store_true", help="Debug mode - verbose")

refdict = {'SPS01_':'Strain 1','SPS02_':'Strain 2','SPS03_':'Strain 3','SPS04_':'Strain 4','SPS05_':'Strain 5','SPS06_':'Strain 6','SPS07_':'Strain 7','SPS08_':'Strain 8',
        'SPS09_':'Strain 9','SPS10_':'Strain 10','SPS11_':'Strain 11','SPS12_':'Strain 12','SPS13_':'Strain 13','SPS14_':'Strain 14','SPS15_':'Strain 15','SPS16_':'Strain 16','SPS17_':'Strain 17',
        'SPS18_':'Strain 18','SPS19_':'Strain 19','SPS20_':'Strain 20','SPS21_':'Strain 21','SPS22_':'Strain 22','SPS23_':'Strain 23','SPS24_':'Strain 24','SPS25_':'Strain 25','SPS26_':'Strain 26',
        'SPS27_':'Strain 27','SPS28_':'Strain 28','SPS29_':'Strain 29','SPS30_':'Strain 30','SP_RS':'TIGR4','SPH_RS':'19A','SPT_RS':'19F','HMPREF0':'TCH84331/19A','MYY_RS':'ST556','SPNA45':'SPNA45',
        'SPP_RS':'P1031','SPNOXC':'OXC141','SPJ_RS':'JJA','SPNINV':'INV200','INV104':'INV104','HMPREF1':'gamPNI0373','SPG_RS':'G54','SPD_RS':'D39','SPCG_R':'CGSP14','SPN23F':'ATCC','SPAP_R':'AP200',
        'SP7058':'70585','SP670_':'670-6B'}


#read fasta file, return dictionary
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles={}
    for i in lines:
        if i[0]==">":
            title=i[:-1]
            titles[title]=""
        else:
            titles[title]+=(i.upper()).replace("\n","")
    return titles

def process_cluster(clusterfile):
    thiscluster=read_sequences(clusterfile)
    return thiscluster

def main():
    opts, args = options.parse_args()
    clusterdir=str(opts.inputdir)
    finaldict = {}
    
    with open(opts.outputfile,"w+") as fout:
        header = ['Product']
        for key,value in refdict.items():
            header.append(value)
        header.append('R6') 
        writer = csv.DictWriter(fout,fieldnames=header)
        writer.writeheader()
        try:
            allClusters=os.listdir(clusterdir)
        except Exception, e:
            print ("Could not find specified clusters directory\n")
            print (str(e))
        
        for cluster in allClusters:
            thiscluster=process_cluster(clusterdir+"/"+cluster)
            thisproduct=thiscluster.keys()[0].split("|")[1]
            forcsv={}
            forcsv['Product'] = thisproduct

            for sequence in thiscluster:
                thisSP = sequence.split("|")[0]
                thisSP = thisSP[1:]

                if 'spr' in thisSP:
                    forcsv['R6'] = thisSP
                elif 'SP_RS' in thisSP:
                    forcsv['TIGR4'] = thisSP
                else:
                    smallSP = thisSP[:6]
                    if smallSP == 'HMPREF':
                        smallSP = thisSP[:7]
                        Strain = refdict[smallSP]
                        forcsv[Strain] = thisSP
                    else:
                        Strain = refdict[smallSP]
                        forcsv[Strain] = thisSP 

            for item in header:
                if item not in forcsv.keys():
                    forcsv[item] = " " 
            writer.writerow(forcsv)

if __name__=="__main__":
    main()
