##usage
#Rscript /jelly_data/yichun/scripts/5formatenrich/05dotplot.R -n Species_tree_gain_loss_dup.txt -a gain_loss_specific.enrichKOG.clean.txt -b gain_loss_specific.enrichKEGG.clean.txt -c gain_loss_specific.enrichGO.clean.txt

#Required the tree structure
##Species_tree_gain_loss_dup.txt

#Required the outputs of 5formatenrich/04gain_loss_specific_enrichment.R
##gain_loss_specific.enrichKOG.clean.txt
##gain_loss_specific.enrichKEGG.clean.txt
##gain_loss_specific.enrichGO.clean.txt

#Output the following

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))

##Create parameters
option_list <- list(
  make_option(c("-n","--nodelabel"), type="character", default=NULL,
              help="Species_tree_gain_loss_dup.txt from gain loss tree' [default %default]",
              dest="speciesorder"),
  make_option(c("-a","--KOG"), type="character", default="gain_loss_specific.enrichKOG.clean.txt",
              help="gain_loss_specific.enrichKOG.txt' [default %default]",
              dest="KOG"),
  make_option(c("-b","--KEGG"), type="character", default="gain_loss_specific.enrichKEGG.clean.txt",
              help="gain_loss_specific.enrichKEGG.txt' [default %default]",
              dest="KEGG"),
  make_option(c("-c","--GO"), type="character", default="gain_loss_specific.enrichGO.clean.txt",
              help="gain_loss_specific.enrichGO.txt' [default %default]",
              dest="GO"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##Manually add files
{
  setwd("D:/centipede/enrich/")
  opt$speciesorder="D:/centipede/cafe/Species_tree_gain_loss_dup.txt"
}

##species order
speciesorder<-read.delim(opt$speciesorder,
                         header = TRUE, stringsAsFactors = FALSE)
speciesorder<-subset(speciesorder, branch.length > 0)
speciesorder.tip<-speciesorder[is.na(as.numeric(gsub("N", "", speciesorder$label)))==TRUE,]
speciesorder.tip<-speciesorder.tip[order(speciesorder.tip$y),]
speciesorder.node<-speciesorder[is.na(as.numeric(gsub("N", "", speciesorder$label)))==FALSE,]
speciesorder.node<-speciesorder.node[order(speciesorder.node$y),]
speciesorder.node<-rbind(speciesorder.tip, speciesorder.node)

####KOG enrich
{
  ##Plot
  plotin<-read.delim(opt$KOG, header = TRUE, sep = "\t")
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotinsep$ratio1=plotinsep$Genenumerator/plotinsep$BGnumerator
  plotinsep$ratio2=plotinsep$Genenumerator/plotinsep$Genedenominator
  plotinsep$ONTOLOGY<-gsub("Information storage and processing","Information storage\nand processing",plotinsep$ONTOLOGY)
  names(plotinsep)[names(plotinsep)=="Groups"]<-"label"
  plotinsep$labelorderBG<-paste("(",plotinsep$BGdenominator,")",plotinsep$label)
  plotinsep$labelordeGE<-paste("(",plotinsep$Genedenominator,")",plotinsep$label)
  
  #Gain
  {
    plotdata<-subset(plotinsep, Type == "Gain")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    a<-ggplot(plotdata, 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "KOG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("KOG.gain.bygroup.tiff", width = 12, height = 10, units = "in", dpi = 300)
    ggsave("KOG.gain.bygroup.png", width = 12, height = 10, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    a<-ggplot(plotdata, 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "KOG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("KOG.gain.bypseudogenome.tiff", width = 12, height = 10, units = "in", dpi = 300)
    ggsave("KOG.gain.bypseudogenome.png", width = 12, height = 10, units = "in", dpi = 300)
  }
  #Loss
  {
    plotdata<-subset(plotinsep, Type == "Loss")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    a<-ggplot(plotdata, 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "KOG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("KOG.loss.bygroup.tiff", width = 12, height = 10, units = "in", dpi = 300)
    ggsave("KOG.loss.bygroup.png", width = 12, height = 10, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    a<-ggplot(plotdata, 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "KOG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("KOG.loss.bypseudogenome.tiff", width = 12, height = 10, units = "in", dpi = 300)
    ggsave("KOG.loss.bypseudogenome.png", width = 12, height = 10, units = "in", dpi = 300)
  }
  #Specific
  {
    plotdata<-subset(plotinsep, Type == "Specific")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.tip, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    a<-ggplot(plotdata, 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "KOG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("KOG.specific.bygroup.tiff", width = 10, height = 10, units = "in", dpi = 300)
    ggsave("KOG.specific.bygroup.png", width = 10, height = 10, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.tip, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    a<-ggplot(plotdata, 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "KOG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("KOG.specific.bygenome.tiff", width = 10, height = 10, units = "in", dpi = 300)
    ggsave("KOG.specific.bygenome.png", width = 10, height = 10, units = "in", dpi = 300)
  }
  
}

####KEGG enrich
{
  ##Plot
  plotin<-read.delim(opt$KEGG, header = TRUE, sep = "\t")
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotinsep$ratio1=plotinsep$Genenumerator/plotinsep$BGnumerator
  plotinsep$ratio2=plotinsep$Genenumerator/plotinsep$Genedenominator
  #plotinsep$ONTOLOGY<-gsub("Information storage and processing","Information storage\nand processing",plotinsep$ONTOLOGY)
  names(plotinsep)[names(plotinsep)=="Groups"]<-"label"
  plotinsep$labelorderBG<-paste("(",plotinsep$BGdenominator,")",plotinsep$label)
  plotinsep$labelordeGE<-paste("(",plotinsep$Genedenominator,")",plotinsep$label)
  
  #Gain
  {
    plotdata<-subset(plotinsep, Type == "Gain")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio2)
    write.table(topgroups, 
                file = "KEGG.gain.bygroup.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "KEGG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("KEGG.gain.bygroup.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("KEGG.gain.bygroup.top10.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio1)
    write.table(topgroups, 
                file = "KEGG.gain.bypseudogenome.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "KEGG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("KEGG.gain.bypseudogenome.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("KEGG.gain.bypseudogenome.top10.png", width = 12, height = 20, units = "in", dpi = 300)
  }
  #Loss
  {
    plotdata<-subset(plotinsep, Type == "Loss")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio2)
    write.table(topgroups, 
                file = "KEGG.loss.bygroup.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "KEGG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("KEGG.loss.bygroup.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("KEGG.loss.bygroup.top10.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio1)
    write.table(topgroups, 
                file = "KEGG.loss.bypseudogenome.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "KEGG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("KEGG.loss.bypseudogenome.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("KEGG.loss.bypseudogenome.top10.png", width = 12, height = 20, units = "in", dpi = 300)
  }
  #Specific
  {
    plotdata<-subset(plotinsep, Type == "Specific")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.tip, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio2)
    write.table(topgroups, 
                file = "KEGG.specific.bygroup.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "KEGG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("KEGG.specific.bygroup.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("KEGG.specific.bygroup.top10.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.tip, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio1)
    write.table(topgroups, 
                file = "KEGG.specific.bygenome.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "KEGG Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      #facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("KEGG.specific.bygenome.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("KEGG.specific.bygenome.top10.png", width = 12, height = 20, units = "in", dpi = 300)
  }
}

####GO enrich
{
  ##Plot
  plotin<-read.delim(opt$GO, header = TRUE, sep = "\t")
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotinsep$ratio1=plotinsep$Genenumerator/plotinsep$BGnumerator
  plotinsep$ratio2=plotinsep$Genenumerator/plotinsep$Genedenominator
  #plotinsep$ONTOLOGY<-gsub("Information storage and processing","Information storage\nand processing",plotinsep$ONTOLOGY)
  names(plotinsep)[names(plotinsep)=="Groups"]<-"label"
  plotinsep$labelorderBG<-paste("(",plotinsep$BGdenominator,")",plotinsep$label)
  plotinsep$labelordeGE<-paste("(",plotinsep$Genedenominator,")",plotinsep$label)
  
  ##Top10
  #Gain
  {
    plotdata<-subset(plotinsep, Type == "Gain")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio2)
    write.table(topgroups, 
                file = "GO.gain.bygroup.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("GO.gain.bygroup.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("GO.gain.bygroup.top10.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio1)
    write.table(topgroups, 
                file = "GO.gain.bypseudogenome.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    topgroups<-subset(topgroups, Count>=5)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("GO.gain.bypseudogenome.top10.tiff", width = 12, height = 40, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GO.gain.bypseudogenome.top10.png", width = 12, height = 40, units = "in", dpi = 300, limitsize = FALSE)
  }
  #Loss
  {
    plotdata<-subset(plotinsep, Type == "Loss")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio2)
    write.table(topgroups, 
                file = "GO.loss.bygroup.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("GO.loss.bygroup.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("GO.loss.bygroup.top10.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio1)
    write.table(topgroups, 
                file = "GO.loss.bypseudogenome.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    topgroups<-subset(topgroups, Count>=5)
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("GO.loss.bypseudogenome.top10.tiff", width = 20, height = 20, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GO.loss.bypseudogenome.top10.png", width = 20, height = 20, units = "in", dpi = 300, limitsize = FALSE)
  }
  #Specific
  {
    plotdata<-subset(plotinsep, Type == "Specific")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.tip, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio2)
    write.table(topgroups, 
                file = "GO.specific.bygroup.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("GO.specific.bygroup.top10.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("GO.specific.bygroup.top10.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.tip, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    topgroups<-plotdata %>% group_by(Cluster) %>% top_n(10, ratio1)
    write.table(topgroups, 
                file = "GO.specific.bygenome.top10.txt",
                row.names = FALSE, sep = "\t", quote = FALSE)
    topgroups<-subset(topgroups, Count>=5)
    a<-ggplot(plotdata[which(plotdata$Description %in% topgroups$Description),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("GO.specific.bygenome.top10.tiff", width = 20, height = 30, units = "in", dpi = 300)
    ggsave("GO.specific.bygenome.top10.png", width = 20, height = 30, units = "in", dpi = 300)
  }
  
  ##Level4
  #Gain
  {
    plotdata<-subset(plotinsep, Type == "Gain")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    level4_terms<-read.delim("D:/3enrichment/GOlevel4.txt", header = TRUE)
    
    plotdata1<-subset(plotdata, goClass %in% level4_terms$goClass)
    topgroups<-plotdata1[which(plotdata1$p.adjust < 2),] %>% group_by(Cluster) %>% top_n(10, ratio2/p.adjust)
    a<-ggplot(subset(plotdata1, p.adjust < 1 & goClass %in% topgroups$goClass), 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("GO.gain.bygroup.L4.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("GO.gain.bygroup.L4.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    plotdata1<-subset(plotdata, goClass %in% level4_terms$goClass)
    topgroups<-plotdata1[which(plotdata1$p.adjust < 2),] %>% group_by(Cluster) %>% top_n(10, ratio2/p.adjust)
    a<-ggplot(subset(plotdata1, p.adjust < 1 & goClass %in% topgroups$goClass), 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("GO.gain.bypseudogenome.L4.tiff", width = 12, height = 20, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GO.gain.bypseudogenome.L4.png", width = 12, height = 20, units = "in", dpi = 300, limitsize = FALSE)
  }
  #Loss
  {
    plotdata<-subset(plotinsep, Type == "Loss")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.node, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    a<-ggplot(plotdata[which(plotdata$goClass %in% level4_terms$goClass),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("GO.loss.bygroup.L4.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("GO.loss.bygroup.L4.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.node, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    a<-ggplot(plotdata[which(plotdata$goClass %in% level4_terms$goClass),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("GO.loss.bypseudogenome.L4.tiff", width = 20, height = 20, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("GO.loss.bypseudogenome.L4.png", width = 20, height = 20, units = "in", dpi = 300, limitsize = FALSE)
  }
  #Specific
  {
    plotdata<-subset(plotinsep, Type == "Specific")
    length(unique(plotdata$label))
    
    labelorderGE<-unique(subset(plotdata, select = c("label","Genedenominator")))
    labelorderGE<-merge(speciesorder.tip, labelorderGE, by = "label", all.x = TRUE, sort = FALSE)
    labelorderGE$labelorderGE<-paste("(",labelorderGE$Genedenominator,")",labelorderGE$label)
    
    a<-ggplot(plotdata[which(plotdata$goClass %in% level4_terms$goClass),], 
              aes(x = labelordeGE, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderGE$labelorderGE, 
                       label = paste(" ",str_replace_all(labelorderGE$labelorderGE,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelordeGE)
    ggsave("GO.specific.bygroup.L4.tiff", width = 12, height = 20, units = "in", dpi = 300)
    ggsave("GO.specific.bygroup.L4.png", width = 12, height = 20, units = "in", dpi = 300)
    
    labelorderBG<-unique(subset(plotdata, select = c("label","BGdenominator")))
    labelorderBG<-merge(speciesorder.tip, labelorderBG, by = "label", all.x = TRUE, sort = FALSE)
    labelorderBG$labelorderBG<-paste("(",labelorderBG$BGdenominator,")",labelorderBG$label)
    
    a<-ggplot(plotdata[which(plotdata$goClass %in% level4_terms$goClass),], 
              aes(x = labelorderBG, y = Description, size = ratio1, colour = p.adjust))+
      labs(title = "GO Enrichment", size = "Ratio", colour = "p.adjust")+
      geom_point(shape = 19)+scale_size_area()+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      scale_x_discrete(limits = labelorderBG$labelorderBG, 
                       label = paste(" ",str_replace_all(labelorderBG$labelorderBG,"_"," ")))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1, face = "italic"), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 0))+
      facet_grid(ONTOLOGY~Type, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 10))
    a
    unique(plotdata$labelorderBG)
    ggsave("GO.specific.bygenome.L4.tiff", width = 20, height = 30, units = "in", dpi = 300)
    ggsave("GO.specific.bygenome.L4.png", width = 20, height = 30, units = "in", dpi = 300)
  }
}