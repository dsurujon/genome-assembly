# genome-assembly
genome assembly for *Streptococcus pneumoniae*

Usage
--------------------
Use this script to assemble complete genomes from a collection of contigs (e.g. obtained from Illumina sequencing)

`assemble_genome.py -i [inputfile] -l [logfile]`

The output will be an unannotated fasta file with the complete genomic sequence. 

The output can then be sent to RAST for annotation.

Genome Statistics
--------------------
Use the `genome_stats.py` script to generate a summary of key features of a genome file. This is usually useful once a newly assembled genome has been annotated

`genome_stats.py -i [inputfile] -g [genome format] -o [outputfile]`

use either `gbk` or `fasta` as the genome format. 

Genome Statistics HTML
-----------------------
This is an updated version of `genome_stats.py` that generates an HTML file that contains more detailed summary statistics including histograms of gene/CDS/RNA feature lengths. 

`genome_stats_html.py -i [inputfile] -g [genome format] -o [output directory] -f [output HTML file]`

The output directory will contain all necessary files for the html file to display properly. 
