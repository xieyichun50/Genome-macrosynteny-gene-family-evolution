##To format the annotation file for single species
##Generate the following files
####Genesorthogrouppair.1v1 -- Gene and orthogroups long table 
####GenesGOpair.1v1 -- Gene and GO long table 
####GenesKOGpair.1v1 -- Gene and KOG long table 
####GenesKEGGpair.1v1 -- Gene and KEGG long table 
####OG.GOpair -- orthogroups and GO long table 
####OG.KOGpair -- orthogroups and KOG long table 
####OG.KEGGpair -- orthogroups and KEGG long table 

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="eggnog_annotation' [default %default]",
              dest="eggnog"),
  make_option(c("-a","--kegg"), type="character", default=NULL,
              help="kegg2name' [default %default]",
              dest="kegg2name"),
  make_option(c("-b","--kog"), type="character", default=NULL,
              help="kog2name' [default %default]",
              dest="kog2name"),
  make_option(c("-c","--go"), type="character", default=NULL,
              help="kog2name' [default %default]",
              dest="kog2name"),
  make_option(c("-o","--output"), type="character", default="out",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

#read in kegg2name
kegg2name <- read.delim(opt$kegg2name,
                        sep = "\t", colClasses = "character")
#read in kog2name
kog2name<-read.delim(opt$kog2name,
                     sep = "\t", colClasses = "character")

#read in GO2name
go2name<-read.delim(opt$go2name,
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"

#read eggnog annotation
eggnog<-read.delim("/home/yichun/huilab/Eurema_hecabe/eggnog/Eurema_hecabe.annotations",
                   header = FALSE)

#KEGG
pathways<-eggnog[,c(1,4)]
pathways<-separate(pathways, V4, c("ko","map"), sep = ",map", remove = FALSE)
pathways<-pathways[,c(1,3)]
names(pathways)[1]="Genes"
names(pathways)[2]="KEGG"
#pathways$KEGG<-gsub("ko", "", pathways$KEGG)
pathways<-subset(pathways, KEGG != "", select = c("Genes", "KEGG"))

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
rm(pairtable, pairtable.new, rcnames, subtable, pathways)
write.table(GenesKEGGpair.1v1, 
            file = "GenesKEGGpair.1v1", 
            row.names = FALSE, sep = "\t")
#KOG
pathways<-eggnog[,c(1,5)]
names(pathways)[1]="Genes"
names(pathways)[2]="KOG"
pathways<-subset(pathways, KOG != "")

#Format GenesKEGGpair
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
rm(pairtable, pairtable.new, rcnames, subtable, pathways)
write.table(GenesKOGpair.1v1, 
            file = "GenesKOGpair.1v1", 
            row.names = FALSE, sep = "\t")

#GO
pathways<-eggnog[,c(1,2)]
rm(eggnog)
names(pathways)[1]="Genes"
names(pathways)[2]="GO"
pathways<-subset(pathways, GO != "")

#Format GenesKEGGpair
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
rm(pairtable, pairtable.new, rcnames, subtable, pathways)
write.table(GenesGOpair.1v1, 
            file = "GenesGOpair.1v1", 
            row.names = FALSE, sep = "\t")

#Orthogroups
orthogroups<-read.delim("/home/yichun/huilab/Eurema_hecabe/longest_prot/OrthoFinder/Results_May05/Orthogroups/Orthogroups.txt",
                        sep = ":", header = FALSE)
names(orthogroups)[1]="orthogroups"
names(orthogroups)[2]="Genes"
#Format Genesorthogrouppair
Genesorthogrouppair.1v1<-matrix(NA, nrow = 1, ncol = 2)
Genesorthogrouppair.1v1<-as.data.frame(Genesorthogrouppair.1v1)

names(Genesorthogrouppair.1v1)[1]="orthogroups"
names(Genesorthogrouppair.1v1)[2]="Genes"

for (i in 1:nrow(orthogroups)) {
  subtable<-orthogroups[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ' ')[[1]]),c(strsplit(subtable$orthogroups[1], ' ')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(orthogroups, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  Genesorthogrouppair.1v1<-rbind(Genesorthogrouppair.1v1, pairtable.new)
}
Genesorthogrouppair.1v1<-subset(Genesorthogrouppair.1v1,
                                is.na(Genesorthogrouppair.1v1$Genes)==FALSE & Genesorthogrouppair.1v1$Genes != "",
                                select = c("orthogroups", "Genes"))
Genesorthogrouppair.1v1<-unique(Genesorthogrouppair.1v1)
rm(pairtable, pairtable.new, rcnames, subtable)
write.table(Genesorthogrouppair.1v1, 
            file = "Genesorthogrouppair.1v1", 
            row.names = FALSE, sep = "\t")


Genesorthogrouppair.1v1<-read.delim("Genesorthogrouppair.1v1", header = TRUE, sep = "\t")
GenesGOpair.1v1<-read.delim("GenesGOpair.1v1", header = TRUE, sep = "\t")
GenesKOGpair.1v1<-read.delim("GenesKOGpair.1v1", header = TRUE, sep = "\t")
GenesKEGGpair.1v1<-read.delim("GenesKEGGpair.1v1", header = TRUE, sep = "\t")

OG.GOpair<-merge(Genesorthogrouppair.1v1, GenesGOpair.1v1, by = "Genes", all.x = TRUE)
OG.GOpair<-subset(OG.GOpair, is.na(GO) == FALSE, select = c("orthogroups", "GO"))
OG.GOpair<-unique(OG.GOpair)
write.table(OG.GOpair, file = "OG.GOpair", row.names = FALSE, sep = "\t")

OG.KOGpair<-merge(Genesorthogrouppair.1v1, GenesKOGpair.1v1, by = "Genes", all.x = TRUE)
OG.KOGpair<-subset(OG.KOGpair, is.na(KOG) == FALSE, select = c("orthogroups", "KOG"))
OG.KOGpair<-unique(OG.KOGpair)
write.table(OG.KOGpair, file = "OG.KOGpair", row.names = FALSE, sep = "\t")

OG.KEGGpair<-merge(Genesorthogrouppair.1v1, GenesKEGGpair.1v1, by = "Genes", all.x = TRUE)
OG.KEGGpair<-subset(OG.KEGGpair, is.na(KEGG) == FALSE, select = c("orthogroups", "KEGG"))
OG.KEGGpair<-unique(OG.KEGGpair)
write.table(OG.KEGGpair, file = "OG.KEGGpair", row.names = FALSE, sep = "\t")

