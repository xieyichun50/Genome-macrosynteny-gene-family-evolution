#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="all-vs-all.blast result to two species' [default %default]",
              dest="input_filename"),
  make_option(c("-o","--output"), type="character", default="ortholog",
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
#opt$input_filename="/home/yichun/huilab/centipede/Orthologues/Orthologues_Trigoniulus_corallinus_hic.gff3.filter.gff3.longest-gene.proteins/Trigoniulus_corallinus_hic.gff3.filter.gff3.longest-gene.proteins__v__Helicorthomorpha_holstii_hic.gff3.filter.gff3.longest-gene.proteins.tsv"
species.pair<-read.delim(opt$input_filename, header = TRUE)

#setwd("/home/yichun/huilab/centipede")
#orthologs<-read.delim("Orthogroups/Orthogroups.tsv", header = TRUE)
#species.pair<-read.delim("Orthologues/Orthologues_Ttu_hic.gff3.filter.gff3.longest-gene.proteins/Ttu_hic.gff3.filter.gff3.longest-gene.proteins__v__Trigoniulus_corallinus_hic.gff3.filter.gff3.longest-gene.proteins.tsv", header = TRUE)
#paralogs<-orthologs[-which(species.pair$Orthogroup %in% orthologs$Orthogroup),]
#paralogs<-subset(paralogs, paralogs$Ttu_hic.gff3.filter.gff3.longest.gene.proteins != "")

names(species.pair)[2]="Species1"
names(species.pair)[3]="Species2"

species.pair.1v1<-matrix(NA, nrow = 1, ncol = 2)
species.pair.1v1<-as.data.frame(species.pair.1v1)
names(species.pair.1v1)[1]="gene1"
names(species.pair.1v1)[2]="gene2"

for (i in 1:nrow(species.pair)) {
  subtable<-species.pair[i,]
  rcnames<-list(c(strsplit(subtable$Species1[1], ', ')[[1]]),c(strsplit(subtable$Species2[1], ', ')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$gene1<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(gene2, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  species.pair.1v1<-rbind(species.pair.1v1, pairtable.new)
}
species.pair.1v1<-subset(species.pair.1v1, is.na(species.pair.1v1$gene1)==FALSE)
write.table(species.pair.1v1, file = paste(opt$output_filename,".",opt$input_filename,".txt", sep = ""), sep = '\t', row.names = FALSE)
