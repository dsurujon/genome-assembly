#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

get_positions<-function(filestr,logstr){
	##read and clean up data
	sink(logstr)
	blast<-read.csv(filestr,header=F)
	cat(paste0("finished reading ",filestr,"\n"))
	names(blast)<-c("contig","cstart","cend","eval","genome","gstart","gend")
	blast$length<-blast$cend-blast$cstart
	blast$group<-interaction(blast$contig, blast$genome)
	blast$forward<-as.integer(blast$gend>blast$gstart)

	ncontig<-length(unique(blast$contig))
	blast<-blast[order(-blast$length),]

	##aggregate by the longest match in each genomeXcontig pair
	blast.longest<-aggregate(.~group, data=blast, FUN=head,1)
	blast.longest$length<-as.numeric(blast.longest$length)

	cat("Summary of longest matches:\n")
	cat(summary(blast.longest$length))
	cat("\n")

	blast.longest$group<-as.character(blast.longest$group)
	
	##find out which genome has the maximum total matches
	blast.longest.totallength<-aggregate(length~genome, data=blast.longest, sum)
	blast.longest.totallength<-blast.longest.totallength[order(-blast.longest.totallength$length),]

	#print(head(blast.longest.totallength))

	#mygenome<-blast.longest.totallength[blast.longest.totallength$length==max(blast.longest.totallength$length),1]
	mygenome<-blast.longest.totallength[1,1]

	#only get the corresponding genome
	blast.longest.gen<-blast.longest[blast.longest$genome==mygenome,]

	j<-1

	while(ncontig!=nrow(blast.longest.gen) & j<21){
		cat(paste0("\nWARNING: NOT ALL CONTIGS COVERED for genome ",as.character(j),"\n"))
		j<-j+1
		mygenome<-blast.longest.totallength[j,1]
		blast.longest.gen<-blast.longest[blast.longest$genome==mygenome,]

	}

	cat("Genome match for longest alignments:\n")
	cat(paste0(blast.longest$group[blast.longest$genome==mygenome][1],"\n"))
	cat("Total match for selected genome:\n")
	cat(paste0(blast.longest.totallength$length[blast.longest.totallength$genome==mygenome],"\n"))
	

	#reorder by genomic position
	blast.longest.gen$min<-pmin(blast.longest.gen$gstart,blast.longest.gen$gend)
	blast.longest.gen<-blast.longest.gen[order(blast.longest.gen$min),c(2,3,4,7,8,10)]
	#contig, cstart, cend, gstart, gend, forward
		

	outfile=paste0(filestr,"_positions.csv")
	write.csv(blast.longest.gen,outfile,quote=F)
	outfile2=paste0(filestr,"_positions.txt")
	write.table(t(blast.longest.gen),outfile2,sep="  ",
			quote=FALSE, row.names=FALSE, col.names=FALSE)

	sink()
}

get_positions(args[1],args[2])

