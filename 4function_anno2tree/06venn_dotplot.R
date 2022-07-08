##Plot venn diagram for groups

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(clusterProfiler))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="orthogroups.tsv' [default %default]",
              dest="orthogroups"),
  make_option(c("-o","--output"), type="character", default="all",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-a","--eggnog"), type="character", default=NULL,
              help="folder with 1v1 functional annotations [default %default]",
              dest="eggnog"),
  make_option(c("-b","--ancestral"), type="character", default=NULL,
              help="folder with ancestral state output [default %default]",
              dest="ancestral"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##Manually add files
{
  setwd("/home/yichun/huilab/longest_prot_without_mite/r8s_filter_lambda1/")
  opt$orthogroups="eggnog/Orthogroups.tsv"
  opt$eggnog="eggnog"
  opt$ancestral="ancestral"
  opt$anno="/home/yichun/3enrichment"
}

orthogroups<-read.delim(opt$orthogroups,
                        sep = "\t", header = TRUE)

ortho.myriapod<-orthogroups

##relax
{
  ortho.millipede<-subset(ortho.myriapod, 
                          Helicorthomorpha_holstii > 0 | 
                            Niponia_nodulosa > 0 |
                            Anaulaciulus_tonginus > 0 |
                            Trigoniulus_corallinus > 0 |
                            Glomeris_maerens > 0)
  ortho.centipede<-subset(ortho.myriapod, 
                          Strigamia_maritima > 0 |
                            Rhysida_immarginata > 0 |
                            Lithobius_niger > 0 |
                            Thereuonema_tuberculata > 0)
  venn.diagram(x = list(Millipede = ortho.millipede$Orthogroup, Centipede = ortho.centipede$Orthogroup),
               filename = "centi2milli/milli2centi.relax.venn.png", 
               imagetype = "png" , units = "in", 
               height = 6, width = 6, resolution = 300,
               cat.pos = 0,
               fill = c("burlywood2", "darkgray"))
}

##Strict
{
  ortho.millipede<-subset(ortho.myriapod, 
                          Helicorthomorpha_holstii > 0 & 
                            Niponia_nodulosa > 0 &
                            Anaulaciulus_tonginus > 0 &
                            Trigoniulus_corallinus > 0 &
                            Glomeris_maerens > 0)
  ortho.centipede<-subset(ortho.myriapod, 
                          Strigamia_maritima > 0 & 
                            Rhysida_immarginata > 0 &
                            Lithobius_niger > 0 &
                            Thereuonema_tuberculata > 0)
  venn.diagram(x = list(Millipede = ortho.millipede$Orthogroup, Centipede = ortho.centipede$Orthogroup),
               filename = "centi2milli/milli2centi.strict.venn.png", 
               imagetype = "png" , units = "in", 
               height = 3, width = 3, resolution = 300,
               cat.pos = 5,
               fill = c("burlywood2", "darkgray"))
}

##pseudo genome
myriapod.pseudogenome<-read.delim(paste0(opt$ancestral,"/pseudo.N9.txt"), header = T)
millipede.pseudogenome<-read.delim(paste0(opt$ancestral,"/pseudo.N14.txt"), header = T)
centipede.pseudogenome<-read.delim(paste0(opt$ancestral,"/pseudo.N13.txt"), header = T)

##get list
ortho.overlap<-ortho.millipede[ortho.millipede$Orthogroup %in% ortho.centipede$Orthogroup, c(1,2)]
names(ortho.overlap)[2]="Groups"
ortho.overlap$Groups= paste0("Common (", nrow(ortho.overlap), ")")
ortho.common<-merge(ortho.overlap, myriapod.pseudogenome, by = "Orthogroup", all = F)

ortho.millipede.spc<-ortho.millipede[-which(ortho.millipede$Orthogroup %in% ortho.overlap$Orthogroup),c(1,2)]
names(ortho.millipede.spc)[2]="Groups"
ortho.millipede.spc$Groups= paste0("Millipede-specific (",nrow(ortho.millipede.spc), ")")
ortho.millipede.spc<-merge(ortho.millipede.spc, millipede.pseudogenome, by = "Orthogroup", all = F)

ortho.centipede.spc<-ortho.centipede[-which(ortho.centipede$Orthogroup %in% ortho.overlap$Orthogroup),c(1,2)]
names(ortho.centipede.spc)[2]="Groups"
ortho.centipede.spc$Groups= paste0("Centipede-specific (",nrow(ortho.centipede.spc), ")")
ortho.centipede.spc<-merge(ortho.centipede.spc, centipede.pseudogenome, by = "Orthogroup", all = F)

rm(ortho.millipede, ortho.centipede, ortho.myriapod)



##function enrich

#read in orthogroup functions
##KEGG
kegg2name<-read.delim(paste(opt$anno, "kegg2name.txt", sep = "/"), 
                      sep = "\t", colClasses = "character")
kegg2name.nh<-read.delim(paste(opt$anno, "kegg2name-nonhuman.txt", sep = "/"), 
                         sep = "\t", colClasses = "character")
kegg2ont<-read.delim(paste(opt$anno,"kegglevel.AC.txt", sep = "/"),
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

myriapod.Orthogroups.KEGG<-read.delim(paste(opt$eggnog, "myriapod.Orthogroups.KEGG.txt", sep = "/"), 
                                      header = TRUE)
myriapod.Orthogroups.KEGG<-merge(myriapod.Orthogroups.KEGG, myriapod.pseudogenome,
                                 by = "Orthogroup", all.x = T)
myriapod.Orthogroups.KEGG<-myriapod.Orthogroups.KEGG[myriapod.Orthogroups.KEGG$KEGG %in% kegg2name.nh$KEGG &
                                                       is.na(myriapod.Orthogroups.KEGG$Genes)==F, c("KEGG","Genes")]
millipede.Orthogroups.KEGG<-read.delim(paste(opt$eggnog, "millipede.Orthogroups.KEGG.txt", sep = "/"), 
                                       header = TRUE)
millipede.Orthogroups.KEGG<-merge(millipede.Orthogroups.KEGG, millipede.pseudogenome,
                                  by = "Orthogroup", all.x = T)
millipede.Orthogroups.KEGG<-millipede.Orthogroups.KEGG[millipede.Orthogroups.KEGG$KEGG %in% kegg2name.nh$KEGG &
                                                       is.na(millipede.Orthogroups.KEGG$Genes)==F, c("KEGG","Genes")]
centipede.Orthogroups.KEGG<-read.delim(paste(opt$eggnog, "centipede.Orthogroups.KEGG.txt", sep = "/"), 
                                       header = TRUE)
centipede.Orthogroups.KEGG<-merge(centipede.Orthogroups.KEGG, centipede.pseudogenome,
                                  by = "Orthogroup", all.x = T)
centipede.Orthogroups.KEGG<-centipede.Orthogroups.KEGG[centipede.Orthogroups.KEGG$KEGG %in% kegg2name.nh$KEGG &
                                                         is.na(centipede.Orthogroups.KEGG$Genes)==F, c("KEGG","Genes")]
##KOG
kog2name<-read.delim(paste(opt$anno, "kog2name.txt", sep = "/"), 
                     sep = "\t", colClasses = "character")
myriapod.Orthogroups.KOG<-read.delim(paste(opt$eggnog, "myriapod.Orthogroups.KOG.txt", sep = "/"), 
                                     header = TRUE)
myriapod.Orthogroups.KOG<-merge(myriapod.Orthogroups.KOG, myriapod.pseudogenome,
                                 by = "Orthogroup", all.x = T)
myriapod.Orthogroups.KOG<-myriapod.Orthogroups.KOG[is.na(myriapod.Orthogroups.KOG$Genes)==F, c("KOG","Genes")]

millipede.Orthogroups.KOG<-read.delim(paste(opt$eggnog, "millipede.Orthogroups.KOG.txt", sep = "/"), 
                                      header = TRUE)
millipede.Orthogroups.KOG<-merge(millipede.Orthogroups.KOG, millipede.pseudogenome,
                                by = "Orthogroup", all.x = T)
millipede.Orthogroups.KOG<-millipede.Orthogroups.KOG[is.na(millipede.Orthogroups.KOG$Genes)==F, c("KOG","Genes")]
centipede.Orthogroups.KOG<-read.delim(paste(opt$eggnog, "centipede.Orthogroups.KOG.txt", sep = "/"), 
                                      header = TRUE)
centipede.Orthogroups.KOG<-merge(centipede.Orthogroups.KOG, centipede.pseudogenome,
                                by = "Orthogroup", all.x = T)
centipede.Orthogroups.KOG<-centipede.Orthogroups.KOG[is.na(centipede.Orthogroups.KOG$Genes)==F, c("KOG","Genes")]

##GO
go2name<-read.delim(paste(opt$anno, "go2name.txt", sep = "/"), 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
myriapod.Orthogroups.GO<-read.delim(paste(opt$eggnog, "myriapod.Orthogroups.GO.txt", sep = "/"), 
                                    header = TRUE)
myriapod.Orthogroups.GO<-merge(myriapod.Orthogroups.GO, myriapod.pseudogenome,
                                by = "Orthogroup", all.x = T)
myriapod.Orthogroups.GO<-myriapod.Orthogroups.GO[is.na(myriapod.Orthogroups.GO$Genes)==F, c("GO","Genes")]
millipede.Orthogroups.GO<-read.delim(paste(opt$eggnog, "millipede.Orthogroups.GO.txt", sep = "/"), 
                                     header = TRUE)
millipede.Orthogroups.GO<-merge(millipede.Orthogroups.GO, millipede.pseudogenome,
                               by = "Orthogroup", all.x = T)
millipede.Orthogroups.GO<-millipede.Orthogroups.GO[is.na(millipede.Orthogroups.GO$Genes)==F, c("GO","Genes")]
centipede.Orthogroups.GO<-read.delim(paste(opt$eggnog, "centipede.Orthogroups.GO.txt", sep = "/"), 
                                     header = TRUE)
centipede.Orthogroups.GO<-merge(centipede.Orthogroups.GO, centipede.pseudogenome,
                               by = "Orthogroup", all.x = T)
centipede.Orthogroups.GO<-centipede.Orthogroups.GO[is.na(centipede.Orthogroups.GO$Genes)==F, c("GO","Genes")]

#GO
{
  GO.myriapod<-compareCluster(Genes ~ Groups, 
                              data = ortho.common, 
                              fun = 'enricher',
                              TERM2GENE = myriapod.Orthogroups.GO,
                              TERM2NAME = go2name,
                              pvalueCutoff = 1,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 1,
                              minGSSize = 1,
                              maxGSSize = 1000000)
  GO.myriapod<-as.data.frame(GO.myriapod)
  GO.millipede<-compareCluster(Genes ~ Groups, 
                               data = ortho.millipede.spc, 
                               fun = 'enricher',
                               TERM2GENE = millipede.Orthogroups.GO,
                               TERM2NAME = go2name,
                               pvalueCutoff = 1,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 10000000)
  GO.millipede<-as.data.frame(GO.millipede)
  GO.centipede<-compareCluster(Genes ~ Groups, 
                               data = ortho.centipede.spc, 
                               fun = 'enricher',
                               TERM2GENE = centipede.Orthogroups.GO,
                               TERM2NAME = go2name,
                               pvalueCutoff = 1,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 10000000)
  GO.centipede<-as.data.frame(GO.centipede)
  
  #Plot
  plotin<-rbind(GO.myriapod, GO.millipede, GO.centipede)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological \n Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular \n Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular \n Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  
  write.table(plotdata, file = "centi2milli/GOenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  level3_terms<-read.delim(paste(opt$anno, "GOlevel3.txt", sep = "/"), header = TRUE)
  level4_terms<-read.delim(paste(opt$anno, "GOlevel4.txt", sep = "/"), header = TRUE)
  ##all GO
  plotdata1<-subset(plotdata, ratio2 >=0.1 & p.adjust < 0.2)
  a<-ggplot(plotdata[which(plotdata$goClass %in% plotdata1$goClass & plotdata$ratio2 > 0.1),], 
            aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Ratio", x = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 1454 )", "Common ( 4880 )", "Centipede-specific ( 1089 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("centi2milli/GO.bygroup.P20.tiff", width = 12, height = 30, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/GO.bygroup.P20.png", width = 12, height = 30, units = "in", dpi = 300, limitsize = FALSE)
  
  ##Level 4  
  
  {
    plotdata1<-subset(plotdata, goClass %in% level4_terms$goClass)
    topgroups<-plotdata1[which(plotdata1$p.adjust < 0.2),] %>% group_by(Cluster) %>% top_n(20, ratio1)
    a<-ggplot(subset(plotdata1, ratio1 > 0 & p.adjust < 1 & goClass %in% topgroups$goClass), 
              aes(x = Groups, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      guides(size = guide_legend(order = 1))+
      scale_x_discrete(limits = c("Millipede-specific (1454)", "Common (4880)", "Centipede-specific (1089)"))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 11))
    a
    ggsave("centi2milli/GO.L4.top20.tiff", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("centi2milli/GO.L4.top20.png", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
    
  }
}

#KEGG
{
  KEGG.myriapod<-compareCluster(Genes ~ Groups, 
                                data = ortho.common, 
                                fun = 'enricher',
                                TERM2GENE = myriapod.Orthogroups.KEGG,
                                TERM2NAME = kegg2name,
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 1,
                                minGSSize = 1,
                                maxGSSize = 200000)
  KEGG.myriapod<-as.data.frame(KEGG.myriapod)
  KEGG.millipede<-compareCluster(Genes ~ Groups, 
                                 data = ortho.millipede.spc, 
                                 fun = 'enricher',
                                 TERM2GENE = millipede.Orthogroups.KEGG,
                                 TERM2NAME = kegg2name,
                                 pvalueCutoff = 1,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 1,
                                 minGSSize = 1,
                                 maxGSSize = 200000)
  KEGG.millipede<-as.data.frame(KEGG.millipede)
  KEGG.centipede<-compareCluster(Genes ~ Groups, 
                                 data = ortho.centipede.spc, 
                                 fun = 'enricher',
                                 TERM2GENE = centipede.Orthogroups.KEGG,
                                 TERM2NAME = kegg2name,
                                 pvalueCutoff = 1,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 1,
                                 minGSSize = 1,
                                 maxGSSize = 200000)
  KEGG.centipede<-as.data.frame(KEGG.centipede)
  
  ##Plot
  plotin<-rbind(KEGG.myriapod, KEGG.millipede, KEGG.centipede)
  
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  

  plotdata<-subset(plotdata, is.na(Description) == FALSE & p.adjust <= 1 & ID != "ko04013")
  plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)
  unique(plotdata$Group)
  write.table(plotdata, file = "centi2milli/KEGG.enrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  plotdata1<-plotdata[plotdata$Count>=2 & is.na(plotdata$ONTOLOGY)==F & plotdata$ONTOLOGY != "Human Diseases",]
  plotdata1$ONTOLOGY<-gsub("Environmental Information Processing", "Environmental Information\nProcessing", plotdata1$ONTOLOGY)
  plotdata1$ONTOLOGY<-gsub("Genetic Information Processing", "Genetic Information\nProcessing", plotdata1$ONTOLOGY)
  
  topgroups<-plotdata1 %>% group_by(Cluster) %>% top_n(20, ratio1*Count)
  a<-ggplot(subset(plotdata1, ratio1 > 0 & ID %in% topgroups$ID ), 
            #  a<-ggplot(subset(plotdata, ratio1 >= 0 & p.adjust < 0.2), 
            aes(x = Groups, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific (1454)", "Common (4880)", "Centipede-specific (1089)"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("centi2milli/KEGG.top20.tiff", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/KEGG.top20.png", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  
  topgroups<-plotdata1[plotdata1$Cluster %in% c("Millipede-specific (1454)","Centipede-specific (1089)") &
                         plotdata1$p.adjust <= 0.2,]
  plotdata1$Groups<-gsub(" ","\n",plotdata1$Groups)
  plotdata1$Groups<-gsub("-","-\n",plotdata1$Groups)
  a<-ggplot(subset(plotdata1, ratio1 > 0 & ID %in% topgroups$ID ), 
            #  a<-ggplot(subset(plotdata, ratio1 >= 0 & p.adjust < 0.2), 
            aes(x = Groups, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-\nspecific\n(1454)", "Common\n(4880)", "Centipede-\nspecific\n(1089)"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("centi2milli/KEGG.specific.tiff", width = 8, height = 8, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/KEGG.specific.png", width = 8, height = 8, units = "in", dpi = 300, limitsize = FALSE)
  
}

#KOG
{
  KOG.myriapod<-compareCluster(Genes ~ Groups, 
                               data = ortho.common, 
                               fun = 'enricher',
                               TERM2GENE = myriapod.Orthogroups.KOG,
                               TERM2NAME = kog2name,
                               pvalueCutoff = 1,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 200000)
  KOG.myriapod<-as.data.frame(KOG.myriapod)
  KOG.millepede<-compareCluster(Genes ~ Groups, 
                                data = ortho.millipede.spc, 
                                fun = 'enricher',
                                TERM2GENE = millipede.Orthogroups.KOG,
                                TERM2NAME = kog2name,
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 1,
                                minGSSize = 1,
                                maxGSSize = 200000)
  KOG.millepede<-as.data.frame(KOG.millepede)
  KOG.centipede<-compareCluster(Genes ~ Groups, 
                                data = ortho.centipede.spc, 
                                fun = 'enricher',
                                TERM2GENE = centipede.Orthogroups.KOG,
                                TERM2NAME = kog2name,
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 1,
                                minGSSize = 1,
                                maxGSSize = 200000)
  KOG.centipede<-as.data.frame(KOG.centipede)
  ##Plot
  plotin<-rbind(KOG.myriapod, KOG.millepede, KOG.centipede)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
  kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  
  write.table(plotdata, file = "centi2milli/KOGenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  a<-ggplot(plotdata, aes(x = Groups, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific (1454)", "Common (4880)", "Centipede-specific (1089)"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("centi2milli/KOG.tiff", width = 12, height = 8, units = "in", dpi = 300)
  ggsave("centi2milli/KOG.png", width = 12, height = 8, units = "in", dpi = 300)
  
}

##Pathview mapping
KEGG2koID<-read.delim(paste(opt$anno,"kegg2koID.txt", sep = "/"),header = TRUE)
millipede.Orthogroups.ko<-read.delim(paste(opt$eggnog, "millipede.Orthogroups.ko.txt", sep = "/"), 
                                     header = TRUE)
millipede.Orthogroups.ko<-millipede.Orthogroups.ko[millipede.Orthogroups.ko$Orthogroup %in% ortho.millipede.spc$Orthogroup,
                                                   c("Orthogroup","ko")]
millipede.Orthogroups.ko<-as.data.frame(xtabs(~ko, millipede.Orthogroups.ko))
names(millipede.Orthogroups.ko)[2]="millipede"

myriapod.Orthogroups.ko<-read.delim(paste(opt$eggnog, "myriapod.Orthogroups.ko.txt", sep = "/"), 
                                     header = TRUE)
myriapod.Orthogroups.ko<-myriapod.Orthogroups.ko[myriapod.Orthogroups.ko$Orthogroup %in% ortho.common$Orthogroup,
                                                   c("Orthogroup","ko")]
myriapod.Orthogroups.ko<-as.data.frame(xtabs(~ko, myriapod.Orthogroups.ko))
names(myriapod.Orthogroups.ko)[2]="myriapod"

pathview.input<-merge(millipede.Orthogroups.ko,myriapod.Orthogroups.ko,
                      by = "ko", all = T)

centipede.Orthogroups.ko<-read.delim(paste(opt$eggnog, "centipede.Orthogroups.ko.txt", sep = "/"), 
                                     header = TRUE)
centipede.Orthogroups.ko<-centipede.Orthogroups.ko[centipede.Orthogroups.ko$Orthogroup %in% ortho.centipede.spc$Orthogroup,
                                                   c("Orthogroup","ko")]
centipede.Orthogroups.ko<-as.data.frame(xtabs(~ko, centipede.Orthogroups.ko))
names(centipede.Orthogroups.ko)[2]="centipede"

pathview.input<-merge(pathview.input,centipede.Orthogroups.ko,
                      by = "ko", all = T)
pathview.input<-pathview.input[is.na(pathview.input$ko)==F,]
pathview.input<-pathview.input[is.na(pathview.input$millipede)==F | is.na(pathview.input$centipede) ==F,]
pathways<-merge(pathview.input, KEGG2koID, by = "ko", all.x = TRUE)
pathways<-as.data.frame(unique(pathways$KEGG))
names(pathways)[1]="KEGG"
pathways<-merge(pathways, kegg2name, by = "KEGG", all.x = TRUE)
names(pathways)[1]="pathways"
pathways<-pathways[is.na(pathways$KEGGname)==F,]
write.table(pathways, "centi2milli/pathways.txt", row.names = F, quote = F, sep = "\t")

row.names(pathview.input)<-pathview.input$ko
pathview.input<-pathview.input[,c("millipede","myriapod","centipede")]
pathview.input[is.na(pathview.input)==T]<-0
write.table(pathview.input, "ko.count.matrix.txt",
            sep = "\t", quote = F, row.names = T)


setwd("centi2milli/")

for (i in 1:nrow(pathways)) {
  pathwayid<-pathways$pathways[i]
  #pathwayid<-"ko04130"
  pv.data<-pathview(gene.data = pathview.input,
                    pathway.id = pathwayid, 
                    species = "ko", 
                    gene.idtype = "KEGG", 
                    limit = list(gene = 1), 
                    bins = list(gene=10), 
                    multi.state = TRUE, 
                    both.dirs = list(gene = F, cpd = T),
                    na.col="transparent", 
                    out.suffix = "milli-centi")
  pv.data<-pv.data$plot.data.gene
}
