
import os
from optparse import OptionParser

options = OptionParser(usage='%prog input output ',
                       description="Specify input gbk file and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file")


def insertline(oldfile,newfile):
	f = open(oldfile)
	lines = f.readlines()
	f.close()
	gcount = 1
	rcount = 1
	locus = "SPS6"       # Must be altered for each genome
	a = open(newfile, 'w+')
	for line in lines[:-1]:
		if "CDS" not in line and line[6:9] != "RNA":
			a.write(line)
		if "CDS" in line:
			a.write(line)
			if gcount < 10:
				a.write("                     /locus_tag=\"" + locus + "_000" + str(gcount) + "\"")
			elif gcount < 100:
				a.write("                     /locus_tag=\"" + locus + "_00" + str(gcount) + "\"")
			elif gcount < 1000:
				a.write("                     /locus_tag=\"" + locus + "_0" + str(gcount) + "\"")
			else:
				a.write("                     /locus_tag=\"" + locus + "_" + str(gcount) + "\"")
			a.write("\n")
			gcount += 1  
		if line[6:9] == "RNA":
			a.write(line)
			if rcount < 10:
				a.write("                     /locus_tag=\"" + locus + "r_00" + str(rcount) + "\"")
			elif rcount < 100:
				a.write("                     /locus_tag=\"" + locus + "r_0" + str(rcount) + "\"")
			else:
				a.write("                     /locus_tag=\"" + locus + "r_" + str(rcount) + "\"")
			a.write("\n")
			rcount += 1  
	a.close()

def main():
	opts, args = options.parse_args()
	insertline(opts.inputfile,opts.outputfile)

if __name__ == '__main__':
	main()     
	
	