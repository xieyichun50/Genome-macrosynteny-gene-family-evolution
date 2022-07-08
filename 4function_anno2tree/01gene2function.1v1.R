##usage
#Rscript /jelly_data/yichun/scripts/5formatenrich/gene2function.1v1.R -i *.eggnog.emapper.annotations -o speciesname
#Output the following
##species.KOG.1v1.txt
##species.KEGG.1v1.txt
##species.GO.1v1.txt
##species.ko.1v1.txt

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="eggnog_annotation' [default %default]",
              dest="eggnog"),
  make_option(c("-o","--output"), type="character", default="out",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

#read in eggnog annotation
#opt$eggnog<-"F:/Myriapod/eggnog/Anaulaciulus_tonginus.eggnog.emapper.annotations"
eggnog<-read.delim(opt$eggnog, header = TRUE, skip = 4)
names(eggnog)[1]="Genes"

#KEGG
pathways<-eggnog[,c("Genes","KEGG_Pathway")]
{
pathways<-separate(pathways, KEGG_Pathway, c("ko","map"), sep = ",map", remove = FALSE)
pathways<-pathways[,c(1,3)]
names(pathways)[1]="Genes"
names(pathways)[2]="KEGG"
pathways<-subset(pathways, KEGG != "-" & KEGG != "" & is.na(KEGG)==FALSE,
                 select = c("Genes", "KEGG"))
pathways<-unique(pathways)
#Format GenesKEGGpair
GenesKEGGpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
GenesKEGGpair.1v1<-as.data.frame(GenesKEGGpair.1v1)
names(GenesKEGGpair.1v1)[1]="Genes"
names(GenesKEGGpair.1v1)[2]="KEGG"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KEGG[1], ',')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(KEGG, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesKEGGpair.1v1<-rbind(GenesKEGGpair.1v1, pairtable.new)
}
GenesKEGGpair.1v1<-subset(GenesKEGGpair.1v1, 
                          is.na(GenesKEGGpair.1v1$Genes)==FALSE, 
                          select = c("KEGG", "Genes"))
GenesKEGGpair.1v1<-unique(GenesKEGGpair.1v1)
}
write.table(GenesKEGGpair.1v1, 
            file = paste(opt$output_filename,".KEGG.1v1.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

##KO
pathways<-eggnog[,c("Genes","KEGG_ko")]
{
names(pathways)[1]="Genes"
  names(pathways)[2]="ko"
  pathways<-subset(pathways, ko != "-" & ko != "" & is.na(ko)==FALSE,
                   select = c("Genes", "ko"))
  pathways<-unique(pathways)
  #Format GenesKEGGpair
  Geneskopair.1v1<-matrix(NA, nrow = 1, ncol = 2)
  Geneskopair.1v1<-as.data.frame(Geneskopair.1v1)
  names(Geneskopair.1v1)[1]="Genes"
  names(Geneskopair.1v1)[2]="ko"
  
  for (i in 1:nrow(pathways)) {
    subtable<-pathways[i,]
    rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$ko[1], ',')[[1]]))
    pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
    pairtable<-as.data.frame(pairtable)
    pairtable$Genes<-rownames(pairtable)
    rownames(pairtable)<-1:nrow(pairtable)
    pairtable<-as.data.frame(pairtable)
    pairtable.new<-pairtable %>% gather(ko, pair, c(1:ncol(pairtable)-1))
    pairtable.new<-pairtable.new[,c(1:2)]
    Geneskopair.1v1<-rbind(Geneskopair.1v1, pairtable.new)
  }
  Geneskopair.1v1<-subset(Geneskopair.1v1, 
                            is.na(Geneskopair.1v1$Genes)==FALSE, 
                            select = c("ko", "Genes"))
  Geneskopair.1v1$ko<-gsub("ko:","",Geneskopair.1v1$ko)
  Geneskopair.1v1<-unique(Geneskopair.1v1)
  
  write.table(Geneskopair.1v1, 
              file = paste(opt$output_filename,".ko.1v1.txt", sep = ""),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  rm(pairtable, pairtable.new, rcnames, subtable, pathways)
}

#KOG
pathways<-eggnog[,c("Genes","COG_category")]
{
names(pathways)[1]="Genes"
names(pathways)[2]="KOG"
pathways<-subset(pathways, KOG != "" & KOG != "-" & is.na(KOG)==FALSE)
pathways<-unique(pathways)

#Format GenesKOGpair
GenesKOGpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
GenesKOGpair.1v1<-as.data.frame(GenesKOGpair.1v1)
names(GenesKOGpair.1v1)[1]="Genes"
names(GenesKOGpair.1v1)[2]="KOG"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KOG[1], '')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(KOG, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesKOGpair.1v1<-rbind(GenesKOGpair.1v1, pairtable.new)
}
GenesKOGpair.1v1<-subset(GenesKOGpair.1v1, 
                         is.na(GenesKOGpair.1v1$Genes)==FALSE, 
                         select = c("KOG", "Genes"))
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1)
}
write.table(GenesKOGpair.1v1, 
            file = paste(opt$output_filename,".KOG.1v1.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

#GO
pathways<-eggnog[,c("Genes", "GOs")]
{
names(pathways)[1]="Genes"
names(pathways)[2]="GO"
pathways<-subset(pathways, GO != "" & GO != "-" & is.na(GO)==FALSE)
pathways<-unique(pathways)

#Format GenesGOGpair
GenesGOpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
GenesGOpair.1v1<-as.data.frame(GenesGOpair.1v1)
names(GenesGOpair.1v1)[1]="Genes"
names(GenesGOpair.1v1)[2]="GO"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$GO[1], ',')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(GO, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesGOpair.1v1<-rbind(GenesGOpair.1v1, pairtable.new)
}
GenesGOpair.1v1<-subset(GenesGOpair.1v1, 
                        is.na(GenesGOpair.1v1$Genes)==FALSE, 
                        select = c("GO", "Genes"))
GenesGOpair.1v1<-unique(GenesGOpair.1v1)
}
write.table(GenesGOpair.1v1, 
            file = paste(opt$output_filename,".GO.1v1.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)
