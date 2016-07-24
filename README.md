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
