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

setwd("D:/centipede")
opt$orthogroups<-"D:/centipede/Orthogroups/Orthogroups.GeneCount.tsv"

orthogroups<-read.delim(opt$orthogroups,
                        sep = "\t", header = TRUE)


ortho.myriapod<-subset(orthogroups, Tachypleus_tridentatus == 0 & Carcinoscorpius_rotundicauda == 0)
ortho.myriapod<-orthogroups
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

##get list
ortho.overlap<-ortho.millipede[ortho.millipede$Orthogroup %in% ortho.centipede$Orthogroup, c(1,13)]
names(ortho.overlap)[2]="Groups"
ortho.overlap$Groups= "Common"
ortho.millipede.spc<-ortho.millipede[-which(ortho.millipede$Orthogroup %in% ortho.overlap$Orthogroup),c(1,13)]
names(ortho.millipede.spc)[2]="Groups"
ortho.millipede.spc$Groups= "Millipede-specific"
ortho.centipede.spc<-ortho.centipede[-which(ortho.centipede$Orthogroup %in% ortho.overlap$Orthogroup),c(1,13)]
names(ortho.centipede.spc)[2]="Groups"
ortho.centipede.spc$Groups= "Centipede-specific"
DEG<-rbind(ortho.overlap,ortho.millipede.spc,ortho.centipede.spc)

##function enrich
##Manually add files
{
  setwd("D:/centipede")
  opt$eggnog="eggnog"
  opt$anno="D:/3enrichment"
}

#read in orthogroup functions
kegg2name<-read.delim(paste(opt$anno, "kegg2name.txt", sep = "/"), 
                      sep = "\t", colClasses = "character")
kegg2name.nh<-read.delim(paste(opt$anno, "kegg2name-nonhuman.txt", sep = "/"), 
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

#annotation for centi+milli
g2.Orthogroups.KOG<-all.Orthogroups.KOG[all.Orthogroups.KOG$Orthogroup %in% DEG$Orthogroup, c(2,1)]
g2.Orthogroups.GO<-all.Orthogroups.GO[all.Orthogroups.GO$Orthogroup %in% DEG$Orthogroup,c(2,1)]
g2.Orthogroups.KEGG<-all.Orthogroups.KEGG[all.Orthogroups.KEGG$Orthogroup %in% DEG$Orthogroup & all.Orthogroups.KEGG$KEGG %in% kegg2name.nh$KEGG, c(2,1)]

#GO
{
  GO.species<-compareCluster(Orthogroup ~ Groups, 
                             data = DEG, 
                             fun = 'enricher',
                             TERM2GENE = g2.Orthogroups.GO,
                             TERM2NAME = go2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 10000000)
  
  #Plot
  plotin<-as.data.frame(GO.species)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[3]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological \n Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular \n Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular \n Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  plotdata$Description<-str_to_title(plotdata$Description)
  plotdata$Description<-gsub(" \\( "," (",plotdata$Description)
  plotdata$Description<-gsub(" \\)",") ",plotdata$Description)
  
  write.table(plotdata, file = "centi2milli/GOenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  level3_terms<-read.delim("D:/3enrichment/GOlevel3.txt", header = TRUE)
  level4_terms<-read.delim("D:/3enrichment/GOlevel4.txt", header = TRUE)
  ##all GO
  plotdata1<-subset(plotdata, ratio2 >=0.1 & p.adjust < 0.2)
  a<-ggplot(plotdata[which(plotdata$goClass %in% plotdata1$goClass & plotdata$ratio2 > 0.1),], 
            aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Ratio", x = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 1378 )", "Common ( 4722 )", "Centipede-specific ( 969 )"))+
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
    topgroups<-plotdata1[which(plotdata1$p.adjust < 2),] %>% group_by(Cluster) %>% top_n(30, ratio1/p.adjust)
    a<-ggplot(subset(plotdata1, ratio1 > 0 & p.adjust < 1 & goClass %in% topgroups$goClass), 
              aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      guides(size = guide_legend(order = 1))+
      scale_x_discrete(limits = c("Millipede-specific ( 1378 )", "Common ( 4722 )", "Centipede-specific ( 969 )"))+
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
    ggsave("centi2milli/GO.bybackground.L4.tiff", width = 12, height = 15, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("centi2milli/GO.bybackground.L4.png", width = 12, height = 15, units = "in", dpi = 300, limitsize = FALSE)
    
    plotdata1<-subset(plotdata, goClass %in% level4_terms$goClass)
    topgroups<-plotdata1[which(plotdata1$p.adjust < 2),] %>% group_by(Cluster) %>% top_n(30, ratio2/p.adjust)
    a<-ggplot(subset(plotdata1, ratio2 > 0 & p.adjust < 2 & goClass %in% topgroups$goClass), 
              aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      guides(size = guide_legend(order = 1))+
      scale_x_discrete(limits = c("Millipede-specific ( 1378 )", "Common ( 4722 )", "Centipede-specific ( 969 )"))+
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
    ggsave("centi2milli/GO.bygroup.L4.tiff", width = 12, height = 15, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("centi2milli/GO.bygroup.L4.png", width = 12, height = 15, units = "in", dpi = 300, limitsize = FALSE)
  }
}

#KEGG
{
  KEGG.species<-compareCluster(Orthogroup ~ Groups, 
                               data = DEG, 
                               fun = 'enricher',
                               TERM2GENE = g2.Orthogroups.KEGG,
                               TERM2NAME = kegg2name,
                               pvalueCutoff = 1,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KEGG.species)
  
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "centi2milli/KEGG.enrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  plotdata<-subset(plotdata, is.na(Description) == FALSE & p.adjust <= 1)
  
  unique(plotdata$Group)
  
  plotdata1<-plotdata
  topgroups<-plotdata1[which(plotdata1$p.adjust < 2),] %>% group_by(Cluster) %>% top_n(20, ratio1/p.adjust)
  a<-ggplot(subset(plotdata1, ratio1 > 0 & p.adjust < 1 & ID %in% topgroups$ID), 
#  a<-ggplot(subset(plotdata, ratio1 >= 0 & p.adjust < 0.2), 
            aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 710 )", "Common ( 2534 )", "Centipede-specific ( 466 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("centi2milli/KEGG.bybackground.top20.tiff", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/KEGG.bybackground.top20.png", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  
  topgroups<-plotdata1[which(plotdata1$p.adjust < 2),] %>% group_by(Cluster) %>% top_n(20, ratio1/p.adjust)
  a<-ggplot(subset(plotdata1, ratio1 > 0 & p.adjust < 1 & ID %in% topgroups$ID), 
  #a<-ggplot(subset(plotdata, ratio2 >= 0 & p.adjust < 0.2),
            aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 710 )", "Common ( 2534 )", "Centipede-specific ( 466 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("centi2milli/KEGG.bygroup.top20.tiff", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/KEGG.bygroup.top20.png", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)

  #top 20
  topgroups<-plotdata %>% group_by(Cluster) %>% top_n(20, ratio1)
  a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
            aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 710 )", "Common ( 2534 )", "Centipede-specific ( 466 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("centi2milli/KEGG.bybackground.top20.r1.tiff", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/KEGG.bybackground.top20.r1.png", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  
  topgroups<-plotdata %>% group_by(Cluster) %>% top_n(20, ratio2)
  a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),],
            aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 710 )", "Common ( 2534 )", "Centipede-specific ( 466 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    #facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("centi2milli/KEGG.bygroup.top20.r2.tiff", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("centi2milli/KEGG.bygroup.top20.r2.png", width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  
}

#KOG
{
  KOG.species<-compareCluster(Orthogroup ~ Groups, 
                              data = DEG, 
                              fun = 'enricher',
                              TERM2GENE = g2.Orthogroups.KOG,
                              TERM2NAME = kog2name,
                              pvalueCutoff = 1,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 1,
                              minGSSize = 1,
                              maxGSSize = 200000)
  ##Plot
  plotin<-as.data.frame(KOG.species)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[3]<-"kogClass"
  kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  plotdata$Description<-paste(plotdata$Description,"(", plotdata$BGnumerator, ")")
  
  write.table(plotdata, file = "centi2milli/KOGenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  a<-ggplot(plotdata, aes(x = Group, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 1546 )", "Common ( 5012 )", "Centipede-specific ( 1203 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("centi2milli/KOG.bybackground.tiff", width = 10, height = 8, units = "in", dpi = 300)
  ggsave("centi2milli/KOG.bybackground.png", width = 10, height = 8, units = "in", dpi = 300)
  
  a<-ggplot(plotdata, aes(x = Group, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    scale_x_discrete(limits = c("Millipede-specific ( 1546 )", "Common ( 5012 )", "Centipede-specific ( 1203 )"))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("centi2milli/KOGp1.bygroup.tiff", width = 10, height = 8, units = "in", dpi = 300)
  ggsave("centi2milli/KOGp1.bygroup.png", width = 10, height = 8, units = "in", dpi = 300)

}
