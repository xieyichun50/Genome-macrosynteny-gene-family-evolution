library(ape)
library(ggplot2)
library(tidytree)
library(phylotools)
library(phangorn)
library(ggtree)
library(flextable)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplotify)
library(eoffice)

setwd("MLtree/")

genelist<-read.delim("genefile", header = FALSE)
for (i in 1:nrow(genelist)) {
  gene<-genelist$V1[i]
  tree<-read.tree(paste0(gene,".align.fa.treefile"))
  
  node.order<-data.frame(node=1:Nnode(tree) + Ntip(tree), node.label = tree$node.label)
  node.order<-separate(node.order, node.label,
                       c("SH-aLRT", "bootstrap"), 
                       sep = "/", convert = FALSE, 
                       remove = FALSE)
  
  node.order$`SH-aLRT`<-as.numeric(node.order$`SH-aLRT`)
  node.order$bootstrap<-as.numeric(node.order$bootstrap)
  
  ##Midpoint
  test<-midpoint(tree)
  mytree<-full_join(test, node.order, by = "node")
  
  a<-ggtree(mytree, size = 1)+
    geom_tiplab(size=6, align=TRUE, linesize=.5)+
    #geom_text2(aes(subset=!isTip,label=bootstrap, hjust=1.5, vjust = -1),size=5, color = "black")+
    geom_nodelab(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.25, nudge_y = 0,
                 size=4, color = "black")+
    geom_treescale()
  
  a
  f=paste0(gene,".rectanglemid.pptx")
  topptx(a,f, width = 6, height = min(node.order$node)*0.4, units = "in")
  
}

