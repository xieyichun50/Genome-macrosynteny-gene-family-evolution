##To format the annotation file for multiple species
##Generate the following files
####Genesorthogrouppair.1v1 -- Gene and orthogroups long table 
####GenesGOpair.1v1 -- Gene and GO long table 
####GenesKOGpair.1v1 -- Gene and KOG long table 
####GenesKEGGpair.1v1 -- Gene and KEGG long table 
####OG.GOpair -- orthogroups and GO long table 
####OG.KOGpair -- orthogroups and KOG long table 
####OG.KEGGpair -- orthogroups and KEGG long table 
setwd("/home/yichun/huilab/Eurema_hecabe/eggnog")
library(dplyr)
library(tidyr)

#read in kegg2name
kegg2name <- read.delim("/home/yichun/tools/Annotation/kegg2name.txt",
                        sep = "\t", colClasses = "character")
#read in kog2name
kog2name<-read.delim("/home/yichun/tools/Annotation/kog2name.txt",
                     sep = "\t", colClasses = "character")

#read in GO2name
go2name<-read.delim("/home/yichun/tools/Annotation/go2name.txt",
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"

specieslist<-read.delim(file = "specieslist", header = FALSE)
specieslist$V1<-as.character(specieslist$V1)

all.GenesGOpair.1v1<-NA
all.GenesKOGpair.1v1<-NA
all.GenesKEGGpair.1v1<-NA
all.Genesorthogrouppair.1v1<-NA

for (j in 1:nrow(specieslist)) {
#read eggnog annotation
  eggnog<-read.delim(paste(specieslist$V1[j], ".eggnog.emapper.annotations", sep = ""), 
                     header = FALSE, skip = 4)
  eggnog<-separate(eggnog, V1, c("Genes", "Transcript"), sep = "-", remove = FALSE)
#############  
#KEGG
pathways<-eggnog[,c(2,12)]
pathways<-separate(pathways, V10, c("ko","map"), sep = ",map", remove = FALSE)
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
            file = paste(specieslist$V1[j], ".GenesKEGGpair.1v1", sep = ""), 
            row.names = FALSE, sep = "\t")
all.GenesKEGGpair.1v1<-rbind(all.GenesKEGGpair.1v1, GenesKEGGpair.1v1)

#############
#KOG
pathways<-eggnog[,c(2,23)]
names(pathways)[1]="Genes"
names(pathways)[2]="KOG"
pathways<-subset(pathways, KOG != "")

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
rm(pairtable, pairtable.new, rcnames, subtable, pathways)
write.table(GenesKOGpair.1v1, 
            file = paste(specieslist$V1[j], ".GenesKOGpair.1v1", sep = ""),
            row.names = FALSE, sep = "\t")
all.GenesKOGpair.1v1<-rbind(all.GenesKOGpair.1v1, GenesKOGpair.1v1)

#############
#GO
pathways<-eggnog[,c(2,9)]
rm(eggnog)
names(pathways)[1]="Genes"
names(pathways)[2]="GO"
pathways<-subset(pathways, GO != "")

#Format GenesGOpair
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
            file = paste(specieslist$V1[j], ".GenesGOpair.1v1", sep = ""),
            row.names = FALSE, sep = "\t")
all.GenesGOpair.1v1<-rbind(all.GenesGOpair.1v1, GenesGOpair.1v1)

#############
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
            file = paste(specieslist$V1[j], ".Genesorthogrouppair.1v1", sep = ""),
            row.names = FALSE, sep = "\t")
all.Genesorthogrouppair.1v1<-rbind(all.Genesorthogrouppair.1v1, Genesorthogrouppair.1v1)
}

##Write to file
all.GenesGOpair.1v1<-subset(all.GenesGOpair.1v1,
                            is.na(Genes)==FALSE)
write.table(all.GenesGOpair.1v1, 
            file = "all.GenesGOpair.1v1",
            row.names = FALSE, sep = "\t")

all.GenesKOGpair.1v1<-subset(all.GenesKOGpair.1v1,
                             is.na(Genes)==FALSE)
write.table(all.GenesKOGpair.1v1, 
            file = "all.GenesKOGpair.1v1",
            row.names = FALSE, sep = "\t")

all.GenesKEGGpair.1v1<-subset(all.GenesKEGGpair.1v1,
                              is.na(Genes)==FALSE)
write.table(all.GenesKEGGpair.1v1, 
            file = "all.GenesKEGGpair.1v1",
            row.names = FALSE, sep = "\t")

all.Genesorthogrouppair.1v1<-subset(all.Genesorthogrouppair.1v1,
                                    is.na(Genes)==FALSE)
write.table(all.Genesorthogrouppair.1v1, 
            file = "all.Genesorthogrouppair.1v1",
            row.names = FALSE, sep = "\t")

OG.GOpair<-merge(all.Genesorthogrouppair.1v1, all.GenesGOpair.1v1, by = "Genes", all.x = TRUE)
OG.GOpair<-subset(OG.GOpair, is.na(GO) == FALSE, select = c("orthogroups", "GO"))
OG.GOpair<-unique(OG.GOpair)
write.table(OG.GOpair, file = "OG.GOpair", row.names = FALSE, sep = "\t")

OG.KOGpair<-merge(all.Genesorthogrouppair.1v1, all.GenesKOGpair.1v1, by = "Genes", all.x = TRUE)
OG.KOGpair<-subset(OG.KOGpair, is.na(KOG) == FALSE, select = c("orthogroups", "KOG"))
OG.KOGpair<-unique(OG.KOGpair)
write.table(OG.KOGpair, file = "OG.KOGpair", row.names = FALSE, sep = "\t")

OG.KEGGpair<-merge(all.Genesorthogrouppair.1v1, all.GenesKEGGpair.1v1, by = "Genes", all.x = TRUE)
OG.KEGGpair<-subset(OG.KEGGpair, is.na(KEGG) == FALSE, select = c("orthogroups", "KEGG"))
OG.KEGGpair<-unique(OG.KEGGpair)
write.table(OG.KEGGpair, file = "OG.KEGGpair", row.names = FALSE, sep = "\t")

