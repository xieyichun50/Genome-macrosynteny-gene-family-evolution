##usage
#Rscript /jelly_data/yichun/scripts/5formatenrich/04gain_loss_specific_enrichment.R -n Species_tree_gain_loss_dup.txt -e eggnog -c cafe -a annotation

#Required the tree structure
##Species_tree_gain_loss_dup.txt

#Required the outputs of 5formatenrich/02orthogroups2function.R
##species.Genesorthogrouppair.1v1.txt
##all.Genesorthogrouppair.1v1.txt
##all.Orthogroups.KOG.txt
##all.Orthogroups.KEGG.txt
##all.Orthogroups.GO.txt 

#Required the outputs of 5formatenrich/03ancestralstate.R
##pseudo.species.txt --pseudogenome of each node/tip of tree
##gain.species.txt --gained genes of each node/tip of tree
##loss.species.txt --lost genes of each node/tip of tree
##genome.species.txt --genomes of species
##specific.species.txt --specific orthologues of species

#Output the following

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(clusterProfiler))

##Create parameters
option_list <- list(
  make_option(c("-n","--nodelabel"), type="character", default=NULL,
              help="Species_tree_gain_loss_dup.txt from gain loss tree' [default %default]",
              dest="speciesorder"),
  make_option(c("-e","--eggnog"), type="character", default=NULL,
              help="folder with 5formatenrich/orthogroups2function.R results' [default %default]",
              dest="eggnog"),
  make_option(c("-c","--cafe"), type="character", default=NULL,
              help="folder with 5formatenrich/03ancestralstate.R results' [default %default]",
              dest="cafe"),
  make_option(c("-a","--anno"), type="character", default=NULL,
              help="folder with annotation headers' [default %default]",
              dest="anno"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="folder to output files' [default %default]",
              dest="output"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##Manually add files
{
  setwd("D:/centipede")
  opt$eggnog="eggnog"
  opt$cafe="cafe"
  opt$speciesorder="Species_tree_gain_loss_dup.txt"
  opt$anno="D:/3enrichment"
  opt$output="enrich"
}

#read in species label
speciesorder<-read.delim(paste(opt$cafe, opt$speciesorder, sep = "/"),
                         header = TRUE, stringsAsFactors = FALSE)

#read in orthogroup functions
kegg2name<-read.delim(paste(opt$anno, "kegg2name.txt", sep = "/"), 
                      sep = "\t", colClasses = "character")
all.Orthogroups.KEGG<-read.delim(paste(opt$eggnog, "all.Orthogroups.KEGG.txt", sep = "/"), 
                                 header = TRUE)

kog2name<-read.delim(paste(opt$anno, "kog2name.txt", sep = "/"), 
                     sep = "\t", colClasses = "character")
all.Orthogroups.KOG<-read.delim(paste(opt$eggnog, "all.Orthogroups.KOG.txt", sep = "/"), 
                                 header = TRUE)

go2name<-read.delim(paste(opt$anno, "go2name.txt", sep = "/"), 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
all.Orthogroups.GO<-read.delim(paste(opt$eggnog, "all.Orthogroups.GO.txt", sep = "/"), 
                                 header = TRUE)

##fast check
{
  for (i in 1:nrow(speciesorder)) {
    if (speciesorder$branch.length[i] <= 0) {
      
    } else {
      speciesname=speciesorder$label[i]
      parentnode<-subset(speciesorder, label == speciesname)
      parentname<-speciesorder$label[which(speciesorder$node==parentnode$parent)]
      #show progress
      cat(i, "speciesname:", speciesname, " ", "parentname:", parentname, "\n")
    }
  }
}

##Create blank matrix
{
  function.all<-matrix(NA, nrow = 1, ncol = 12, 
                  dimnames = list(NA, c("Cluster","Groups",
                                        "ID","Description",
                                        "GeneRatio","BgRatio",
                                        "pvalue","p.adjust",
                                        "qvalue","geneID",
                                        "Count","Type")))
  function.all<-as.data.frame(function.all)
  KOG.all<-function.all
  GO.all<-function.all
  KEGG.all<-function.all
}

##Statistics
{
  for (i in 12:nrow(speciesorder)) {
    if (speciesorder$branch.length[i] <= 0) {
      
    } else {
      speciesname=speciesorder$label[i]
      parentnode<-subset(speciesorder, label == speciesname)
      parentname<-speciesorder$label[which(speciesorder$node==parentnode$parent)]
      #show progress
      cat(i, "speciesname:", speciesname, " ", "parentname:", parentname, "\n")
      
      ##Enrich gained genes
      {
        gaingenes<-read.delim(paste(opt$cafe, "/gain.",speciesname,".txt", sep = ""), 
                              header = TRUE)
        gaingenes$Groups=speciesname
        
        currentgenome<-read.delim(paste(opt$cafe, "/", "pseudo.", speciesname, ".txt", sep = ""), 
                                  header = TRUE)
        #KOG
        {
          currentgenome.KOG<-merge(currentgenome, all.Orthogroups.KOG, by = "Orthogroup", all.x = TRUE)
          currentgenome.KOG<-subset(currentgenome.KOG, is.na(KOG)==FALSE, select = c("KOG", "Genes"))
          
          KOG.species<-compareCluster(Genes ~ Groups, 
                                      data = gaingenes, 
                                      fun = 'enricher',
                                      TERM2GENE = currentgenome.KOG,
                                      TERM2NAME = kog2name,
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 1,
                                      maxGSSize = 200000)
          KOG.species<-as.data.frame(KOG.species)
          KOG.species$Type="Gain"
          KOG.all<-rbind(KOG.all, KOG.species)
          rm(currentgenome.KOG, KOG.species)
        }
        #KEGG
        {
          currentgenome.KEGG<-merge(currentgenome, all.Orthogroups.KEGG, by = "Orthogroup", all.x = TRUE)
          currentgenome.KEGG<-subset(currentgenome.KEGG, is.na(KEGG)==FALSE, select = c("KEGG", "Genes"))
          
          KEGG.species<-compareCluster(Genes ~ Groups, 
                                       data = gaingenes, 
                                       fun = 'enricher',
                                       TERM2GENE = currentgenome.KEGG,
                                       TERM2NAME = kegg2name,
                                       pvalueCutoff = 1,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 1,
                                       minGSSize = 1,
                                       maxGSSize = 200000)
          KEGG.species<-as.data.frame(KEGG.species)
          KEGG.species$Type="Gain"
          KEGG.all<-rbind(KEGG.all, KEGG.species)
          rm(currentgenome.KEGG, KEGG.species)
        }
        #GO
        {
          currentgenome.GO<-merge(currentgenome, all.Orthogroups.GO, by = "Orthogroup", all.x = TRUE)
          currentgenome.GO<-subset(currentgenome.GO, is.na(GO)==FALSE, select = c("GO", "Genes"))
          
          GO.species<-compareCluster(Genes ~ Groups, 
                                     data = gaingenes, 
                                     fun = 'enricher',
                                     TERM2GENE = currentgenome.GO,
                                     TERM2NAME = go2name,
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 1,
                                     minGSSize = 1,
                                     maxGSSize = 10000000)
          GO.species<-as.data.frame(GO.species)
          GO.species$Type="Gain"
          GO.all<-rbind(GO.all, GO.species)
          rm(currentgenome.GO, GO.species)
        }
        
        rm(gaingenes, currentgenome)
      }
      
      ##Enrich lost genes
      {
        lossgenes<-read.delim(paste(opt$cafe, "/loss.",speciesname,".txt", sep = ""), 
                              header = TRUE)
        lossgenes$Groups=speciesname
        
        parentnode<-subset(speciesorder, label == speciesname)
        parentname<-speciesorder$label[which(speciesorder$node==parentnode$parent)]
        currentgenome<-read.delim(paste(opt$cafe, "/", "pseudo.", parentname, ".txt", sep = ""), 
                                  header = TRUE)
        #KOG
        {
          currentgenome.KOG<-merge(currentgenome, all.Orthogroups.KOG, by = "Orthogroup", all.x = TRUE)
          currentgenome.KOG<-subset(currentgenome.KOG, is.na(KOG)==FALSE, select = c("KOG", "Genes"))
          
          KOG.species<-compareCluster(Genes ~ Groups, 
                                      data = lossgenes, 
                                      fun = 'enricher',
                                      TERM2GENE = currentgenome.KOG,
                                      TERM2NAME = kog2name,
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 1,
                                      maxGSSize = 200000)
          KOG.species<-as.data.frame(KOG.species)
          KOG.species$Type="Loss"
          KOG.all<-rbind(KOG.all, KOG.species)
          rm(currentgenome.KOG, KOG.species)
        }
        #KEGG
        {
          currentgenome.KEGG<-merge(currentgenome, all.Orthogroups.KEGG, by = "Orthogroup", all.x = TRUE)
          currentgenome.KEGG<-subset(currentgenome.KEGG, is.na(KEGG)==FALSE, select = c("KEGG", "Genes"))
          
          KEGG.species<-compareCluster(Genes ~ Groups, 
                                       data = lossgenes, 
                                       fun = 'enricher',
                                       TERM2GENE = currentgenome.KEGG,
                                       TERM2NAME = kegg2name,
                                       pvalueCutoff = 1,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 1,
                                       minGSSize = 1,
                                       maxGSSize = 200000)
          KEGG.species<-as.data.frame(KEGG.species)
          KEGG.species$Type="Loss"
          KEGG.all<-rbind(KEGG.all, KEGG.species)
          rm(currentgenome.KEGG, KEGG.species)
        }
        #GO
        {
          currentgenome.GO<-merge(currentgenome, all.Orthogroups.GO, by = "Orthogroup", all.x = TRUE)
          currentgenome.GO<-subset(currentgenome.GO, is.na(GO)==FALSE, select = c("GO", "Genes"))
          
          GO.species<-compareCluster(Genes ~ Groups, 
                                     data = lossgenes, 
                                     fun = 'enricher',
                                     TERM2GENE = currentgenome.GO,
                                     TERM2NAME = go2name,
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 1,
                                     minGSSize = 1,
                                     maxGSSize = 10000000)
          GO.species<-as.data.frame(GO.species)
          GO.species$Type="Loss"
          GO.all<-rbind(GO.all, GO.species)
          rm(currentgenome.GO, GO.species)
        }
        
        rm(lossgenes, currentgenome, parentnode, parentname)
      }
      
      ##Enrich species-specific genes
      if (is.na(as.numeric(gsub("N", "", speciesname)))){
        spcfgenes<-read.delim(paste(opt$cafe, "/specific.",speciesname,".txt", sep = ""), 
                              header = TRUE)
        spcfgenes$Groups=speciesname
        
        currentgenome<-read.delim(paste(opt$cafe, "/", "genome.", speciesname, ".txt", sep = ""), 
                                  header = TRUE)
        currentgenome.KOG<-merge(spcfgenes, all.Orthogroups.KOG, by = "Orthogroup", all.x = TRUE)
        
        #KOG
        {
          currentgenome.KOG<-merge(currentgenome, all.Orthogroups.KOG, by = "Orthogroup", all.x = TRUE)
          currentgenome.KOG<-subset(currentgenome.KOG, is.na(KOG)==FALSE, select = c("KOG", "Genes"))
          
          KOG.species<-compareCluster(Genes ~ Groups, 
                                      data = spcfgenes, 
                                      fun = 'enricher',
                                      TERM2GENE = currentgenome.KOG,
                                      TERM2NAME = kog2name,
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 1,
                                      minGSSize = 1,
                                      maxGSSize = 200000)
          KOG.species<-as.data.frame(KOG.species)
          KOG.species$Type="Specific"
          KOG.all<-rbind(KOG.all, KOG.species)
          rm(currentgenome.KOG, KOG.species)
        }
        #KEGG
        {
          currentgenome.KEGG<-merge(currentgenome, all.Orthogroups.KEGG, by = "Orthogroup", all.x = TRUE)
          currentgenome.KEGG<-subset(currentgenome.KEGG, is.na(KEGG)==FALSE, select = c("KEGG", "Genes"))
          
          KEGG.species<-compareCluster(Genes ~ Groups, 
                                       data = spcfgenes, 
                                       fun = 'enricher',
                                       TERM2GENE = currentgenome.KEGG,
                                       TERM2NAME = kegg2name,
                                       pvalueCutoff = 1,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 1,
                                       minGSSize = 1,
                                       maxGSSize = 200000)
          KEGG.species<-as.data.frame(KEGG.species)
          KEGG.species$Type="Specific"
          KEGG.all<-rbind(KEGG.all, KEGG.species)
          rm(currentgenome.KEGG, KEGG.species)
        }
        #GO
        {
          currentgenome.GO<-merge(currentgenome, all.Orthogroups.GO, by = "Orthogroup", all.x = TRUE)
          currentgenome.GO<-subset(currentgenome.GO, is.na(GO)==FALSE, select = c("GO", "Genes"))
          
          GO.species<-compareCluster(Genes ~ Groups, 
                                     data = spcfgenes, 
                                     fun = 'enricher',
                                     TERM2GENE = currentgenome.GO,
                                     TERM2NAME = go2name,
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 1,
                                     minGSSize = 1,
                                     maxGSSize = 10000000)
          GO.species<-as.data.frame(GO.species)
          GO.species$Type="Specific"
          GO.all<-rbind(GO.all, GO.species)
          rm(currentgenome.GO, GO.species)
        }
        
        rm(spcfgenes, currentgenome)
      } else {
        cat("skip specific gene analysis of ", speciesname)
      }
    }
  }
}

##Write to file
KOG.all<-KOG.all[which(is.na(KOG.all$Description)==FALSE),]
names(KOG.all)[3]<-"kogClass"
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
KOG.all<-merge(KOG.all, kog2name, by = c("kogClass"), all.x = TRUE)
KOG.all<-unique(KOG.all)

write.table(KOG.all, 
            file = paste(opt$output,"/gain_loss_specific.enrichKOG.txt", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

KOG.all.clean<-KOG.all[,-which(names(KOG.all)=="geneID")]
write.table(KOG.all.clean, 
            file = paste(opt$output,"/gain_loss_specific.enrichKOG.clean.txt", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

GO.all<-GO.all[which(is.na(GO.all$Description)==FALSE),]
colnames(GO.all)[3]<-"goClass"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

GO.all<-merge(GO.all, go2name, by = c("goClass"), all.x = TRUE)
GO.all<-unique(GO.all)

write.table(GO.all, 
            file = paste(opt$output,"/gain_loss_specific.enrichGO.txt", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

GO.all.clean<-GO.all[,-which(names(GO.all)=="geneID")]
write.table(GO.all.clean, 
            file = paste(opt$output,"/gain_loss_specific.enrichGO.clean.txt", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

KEGG.all<-KEGG.all[which(is.na(KEGG.all$Description)==FALSE),]
KEGG.all<-unique(KEGG.all)

write.table(KEGG.all, 
            file = paste(opt$output,"/gain_loss_specific.enrichKEGG.txt", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

KEGG.all.clean<-KEGG.all[,-which(names(KEGG.all)=="geneID")]
write.table(KEGG.all.clean, 
            file = paste(opt$output,"/gain_loss_specific.enrichKEGG.clean.txt", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)
