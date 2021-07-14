#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="species.longest-gene.gtf [default %default]",
              dest="input_filename"),
  make_option(c("-o","--output"), type="character", default="formatted",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-a", "--xspecies-name"), type="character", default="x.species",
              help="Reference name show at xlab [default %default]",
              dest="name_x"),
  make_option(c("-b", "--yspecies-name"), type="character", default="y.species",
              help="query name show at ylab [default %default]",
              dest="name_y")
)
options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

####Format genome gtf
gtf<-read.table(opt$input_filename, stringsAsFactors = F, skip = 29)
gtf<-gtf[,c(1,4,5,10)]
names(gtf)[1]="chr"
names(gtf)[2]="start"
names(gtf)[3]="end"
names(gtf)[4]="gene"
genelist<-unique(gtf$gene)
gff3<-gtf[1,]
gff3[,1:4]<-NA
for (i in 1:length(genelist)) {
  x.table<-subset(gtf, gene == genelist[i])
  row.names(x.table)=1:nrow(x.table)
  x.table[1,2]<-min(x.table$start)
  x.table[1,3]<-max(x.table$end)
  x.table<-x.table[1,]
  gff3<-rbind(gff3, x.table)
}
gff3<-subset(gff3, is.na(chr)==FALSE)

####for (x,y) position
gff3<-gff3[,c("chr", "start", "end", "gene")]
write.table(gff3, file = paste(opt$output_filename,".",opt$input_filename, sep = ""), sep = '\t', row.names = FALSE)
