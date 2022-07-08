##usage
#Rscript /jelly_data/yichun/scripts/5formatenrich/orthogroups2function.R -i OrthoFinder/*/Orthogroups.tsv 

#Required the outputs of 5formatenrich/gene2function.1v1.R
##species.KOG.1v1.txt
##species.KEGG.1v1.txt
##species.GO.1v1.txt
##species.ko.1v1.txt

#Output the following
##myriapod.Genesorthogrouppair.1v1.txt
##myriapod.Orthogroups.KOG.txt
##myriapod.Orthogroups.KEGG.txt
##myriapod.Orthogroups.GO.txt 
##myriapod.Orthogroups.ko.txt

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="orthogroups' [default %default]",
              dest="orthogroups"),
  make_option(c("-o","--output"), type="character", default="all",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

#read in orthogroups
#setwd("D:/centipede")
#opt$orthogroups<-"D:/centipede/Orthogroups/Orthogroups.tsv"

orthogroups<-read.delim(opt$orthogroups,
                        sep = "\t", header = TRUE)
##Unmask and Run One of the following
{
  #ref="myriapod"
  #specieslist<-c("Orthogroup","Helicorthomorpha_holstii",
  #               "Niponia_nodulosa","Anaulaciulus_tonginus",
  #               "Trigoniulus_corallinus","Glomeris_maerens",
  #               "Strigamia_maritima","Rhysida_immarginata",
  #               "Lithobius_niger","Thereuonema_tuberculata")
}

{
  #ref="millipede"
  #specieslist<-c("Orthogroup","Helicorthomorpha_holstii",
  #               "Niponia_nodulosa","Anaulaciulus_tonginus",
  #               "Trigoniulus_corallinus","Glomeris_maerens")
}

{
  #ref="centipede"
  #specieslist<-c("Orthogroup",
  #               "Strigamia_maritima","Rhysida_immarginata",
  #               "Lithobius_niger","Thereuonema_tuberculata")
}

orthogroups<-orthogroups[,names(orthogroups) %in% specieslist]
##Format Genesorthogrouppair
{
  ##all.Genesorthogrouppair.1v1
  all.Genesorthogrouppair.1v1<-matrix(NA, nrow = 0, ncol = 3)
  all.Genesorthogrouppair.1v1<-as.data.frame(all.Genesorthogrouppair.1v1)
  
  names(all.Genesorthogrouppair.1v1)[1]="Orthogroup"
  names(all.Genesorthogrouppair.1v1)[2]="Genes"
  names(all.Genesorthogrouppair.1v1)[3]="Species"
  
  ##all.Orthogroups.KEGG
  all.Orthogroups.KEGG<-matrix(NA, nrow = 0, ncol = 2)
  all.Orthogroups.KEGG<-as.data.frame(all.Orthogroups.KEGG)
  names(all.Orthogroups.KEGG)[1]="Orthogroup"
  names(all.Orthogroups.KEGG)[2]="KEGG"
  
  ##all.Orthogroups.KOG
  all.Orthogroups.KOG<-all.Orthogroups.KEGG
  names(all.Orthogroups.KOG)[2]="KOG"
  
  ##all.Orthogroups.GO
  all.Orthogroups.GO<-all.Orthogroups.KEGG
  names(all.Orthogroups.GO)[2]="GO"
  
  ##all.Orthogroups.ko
  all.Orthogroups.ko<-all.Orthogroups.KEGG
  names(all.Orthogroups.ko)[2]="ko"
}

##Generate orthogroup-gene-function pair tables
for (j in 2:ncol(orthogroups)) {
  subspecies<-orthogroups[,c(1,j)]
  speciesname=names(subspecies)[2]
  cat(paste0(speciesname,"\n"))
  names(subspecies)[2]="Genes"
  subspecies<-subset(subspecies, is.na(Genes)==FALSE & Genes != "")
  
  #Genes orthogroup pair
  {
    Genesorthogrouppair.1v1<-matrix(NA, nrow = 0, ncol = 2)
    Genesorthogrouppair.1v1<-as.data.frame(Genesorthogrouppair.1v1)
    
    names(Genesorthogrouppair.1v1)[1]="Orthogroup"
    names(Genesorthogrouppair.1v1)[2]="Genes"
    
    for (i in 1:nrow(subspecies)) {
      subtable<-subspecies[i,]
      rcnames<-list(c(strsplit(subtable$Genes[1], ', ')[[1]]),c(strsplit(subtable$Orthogroup[1], ' ')[[1]]))
      pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
      pairtable<-as.data.frame(pairtable)
      pairtable$Genes<-rownames(pairtable)
      rownames(pairtable)<-1:nrow(pairtable)
      pairtable<-as.data.frame(pairtable)
      pairtable.new<-pairtable %>% gather(Orthogroup, pair, c(1:ncol(pairtable)-1))
      pairtable.new<-pairtable.new[,c(1:2)]
      Genesorthogrouppair.1v1<-rbind(Genesorthogrouppair.1v1, pairtable.new)
    }
    
    Genesorthogrouppair.1v1<-subset(Genesorthogrouppair.1v1,
                                    is.na(Genesorthogrouppair.1v1$Genes)==FALSE & Genesorthogrouppair.1v1$Genes != "",
                                    select = c("Orthogroup", "Genes"))
    Genesorthogrouppair.1v1<-unique(Genesorthogrouppair.1v1)
    Genesorthogrouppair.1v1$Species<-speciesname
    rm(pairtable, pairtable.new, rcnames, subtable)
    write.table(Genesorthogrouppair.1v1, 
                file = paste(speciesname, ".Genesorthogrouppair.1v1.txt", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
    all.Genesorthogrouppair.1v1<-rbind(all.Genesorthogrouppair.1v1, Genesorthogrouppair.1v1)
  }
  
  #KEGG
  {
    GenesKEGGpair.1v1<-read.delim(file = paste(speciesname, ".fa.KEGG.1v1.txt", sep = ""),
                                  header = TRUE, sep = "\t")
    Orthogroups.KEGG<-merge(Genesorthogrouppair.1v1, GenesKEGGpair.1v1, 
                            by = "Genes", all.x = TRUE)
    Orthogroups.KEGG<-subset(Orthogroups.KEGG, is.na(KEGG)==FALSE, select = c("Orthogroup","KEGG"))
    Orthogroups.KEGG<-unique(Orthogroups.KEGG[order(Orthogroups.KEGG$Orthogroup),])
    all.Orthogroups.KEGG<-rbind(all.Orthogroups.KEGG, Orthogroups.KEGG)
    all.Orthogroups.KEGG<-unique(all.Orthogroups.KEGG)
    rm(GenesKEGGpair.1v1, Orthogroups.KEGG)
  }
  
  #KOG
  {
    GenesKOGpair.1v1<-read.delim(file = paste(speciesname, ".fa.KOG.1v1.txt", sep = ""),
                                 header = TRUE, sep = "\t")
    Orthogroups.KOG<-merge(Genesorthogrouppair.1v1, GenesKOGpair.1v1, 
                           by = "Genes", all.x = TRUE)
    Orthogroups.KOG<-subset(Orthogroups.KOG, is.na(KOG)==FALSE, select = c("Orthogroup","KOG"))
    Orthogroups.KOG<-unique(Orthogroups.KOG[order(Orthogroups.KOG$Orthogroup),])
    all.Orthogroups.KOG<-rbind(all.Orthogroups.KOG, Orthogroups.KOG)
    all.Orthogroups.KOG<-unique(all.Orthogroups.KOG)
    rm(GenesKOGpair.1v1, Orthogroups.KOG)
  }
  
  #GO
  {
    GenesGOpair.1v1<-read.delim(file = paste(speciesname, ".fa.GO.1v1.txt", sep = ""),
                                header = TRUE, sep = "\t")
    Orthogroups.GO<-merge(Genesorthogrouppair.1v1, GenesGOpair.1v1, 
                          by = "Genes", all.x = TRUE)
    Orthogroups.GO<-subset(Orthogroups.GO, is.na(GO)==FALSE, select = c("Orthogroup","GO"))
    Orthogroups.GO<-unique(Orthogroups.GO[order(Orthogroups.GO$Orthogroup),])
    all.Orthogroups.GO<-rbind(all.Orthogroups.GO, Orthogroups.GO)
    all.Orthogroups.GO<-unique(all.Orthogroups.GO)
    rm(GenesGOpair.1v1, Orthogroups.GO)
  }
  
  #ko
  {
    Geneskopair.1v1<-read.delim(file = paste(speciesname, ".fa.ko.1v1.txt", sep = ""),
                                header = TRUE, sep = "\t")
    Orthogroups.ko<-merge(Genesorthogrouppair.1v1, Geneskopair.1v1, 
                          by = "Genes", all.x = TRUE)
    Orthogroups.ko<-subset(Orthogroups.ko, is.na(ko)==FALSE, select = c("Orthogroup","ko"))
    Orthogroups.ko<-unique(Orthogroups.ko[order(Orthogroups.ko$Orthogroup),])
    all.Orthogroups.ko<-rbind(all.Orthogroups.ko, Orthogroups.ko)
    all.Orthogroups.ko<-unique(all.Orthogroups.ko)
    rm(Geneskopair.1v1, Orthogroups.ko)
  }
}

##Write to file
all.Genesorthogrouppair.1v1<-subset(all.Genesorthogrouppair.1v1,
                                    is.na(Genes)==FALSE)
all.Genesorthogrouppair.1v1<-unique(all.Genesorthogrouppair.1v1)
write.table(all.Genesorthogrouppair.1v1, 
            file = paste0(ref,".Genesorthogrouppair.1v1.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)

all.Orthogroups.KEGG<-subset(all.Orthogroups.KEGG,
                             is.na(KEGG)==FALSE)
all.Orthogroups.KEGG<-unique(all.Orthogroups.KEGG)
write.table(all.Orthogroups.KEGG, 
            file = paste0(ref,".Orthogroups.KEGG.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)

all.Orthogroups.KOG<-subset(all.Orthogroups.KOG,
                            is.na(KOG)==FALSE)
all.Genesorthogrouppair.1v1<-unique(all.Genesorthogrouppair.1v1)
write.table(all.Orthogroups.KOG, 
            file = paste0(ref,".Orthogroups.KOG.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)

all.Orthogroups.GO<-subset(all.Orthogroups.GO,
                           is.na(GO)==FALSE)
all.Orthogroups.GO<-unique(all.Orthogroups.GO)
write.table(all.Orthogroups.GO, 
            file = paste0(ref,".Orthogroups.GO.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)

all.Orthogroups.ko<-subset(all.Orthogroups.ko,
                           is.na(ko)==FALSE)
all.Orthogroups.ko<-unique(all.Orthogroups.ko)
write.table(all.Orthogroups.ko, 
            file = paste0(ref,".Orthogroups.ko.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)
