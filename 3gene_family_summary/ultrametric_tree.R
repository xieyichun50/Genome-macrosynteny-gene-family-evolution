library(phytools)
library(phangorn) 
library(ape)
library(ggplot2)
library(tidytree)
library(ggtree)
library(flextable)
library(tidyr)
library(dplyr)
library(stringr)
library(svglite)
library(ggplotify)
library(eoffice)
setwd("D:/coral/")
tree_orthofinder = read.newick(file = "SpeciesTree_rooted_node_labels.txt",text)

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls"){
    tree<-nnls.tree(cophenetic(tree),tree,
                    rooted=TRUE,trace=0)
  } else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

tree_orthofinder_ultra = force.ultrametric(tree_orthofinder)
#tree_orthofinder_ultra = force.ultrametric(tree_orthofinder, method = "extend")

is.ultrametric(tree_orthofinder_ultra)
write.tree(tree_orthofinder_ultra, file = "tree_ultrametric.txt", append = FALSE, digits = 10, tree.names = FALSE)

#test edge lengths
tree_orthofinder_ultra<-reorder(tree_orthofinder_ultra)
h.nnls<-rowMeans(nodeHeights(tree_orthofinder_ultra))
plot(h.nnls,tree_orthofinder$edge.length-tree_orthofinder_ultra$edge.length,pch=21,
     #ylim=c(-1e-6,1e-6),
     bg="grey",cex=1.5,xlab="edge height",
     ylab="difference between input & output edge lengths",
     main="force.ultrametric(...,method=\"nnls\")")

sum((tree_orthofinder$edge.length-tree_orthofinder_ultra$edge.length)^2)


p<-ggtree(tree_orthofinder, size = 1)+
  geom_tiplab(aes(label=paste("     ",str_replace_all(label,"_"," ")), sep = ""),
              offset = 0.05,
              font = "italic",
              size = 6, color = "black")+
  geom_nodelab(aes(label=label),
             size=6, color = "black")
  xlim(0, 2)

p
ggsave("SpeciesTree_rooted_node_labels.png", width = 20, height = 10, units = "in", dpi = 300)
ggsave("SpeciesTree_rooted_node_labels.tiff", width = 10, height = 10, units = "in", dpi = 300)


p<-ggtree(tree_orthofinder_ultra, size = 1)+
  geom_tiplab(aes(label=paste("     ",str_replace_all(label,"_"," ")), sep = ""),
              offset = 0.05,
              font = "italic",
              size = 6, color = "black")+
  geom_nodelab(aes(label=label),
               size=6, color = "black")
  xlim(0, 2)

p
ggsave("tree_ultrametric.png", width = 20, height = 10, units = "in", dpi = 300)
ggsave("tree_ultrametric.tiff", width = 10, height = 10, units = "in", dpi = 300)
