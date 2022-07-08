library(ape)
library(ggplot2)
library(tidytree)
library(ggtree)
library(flextable)
library(tidyr)
library(dplyr)
library(stringr)
#library(svglite)
library(ggplotify)
library(eoffice)

setwd("/home/yichun/huilab/")
setwd("longest_prot_without_human/r8s_filter_lambda1/")
setwd("longest_prot_without_mite/r8s_filter_lambda1/")

##cafe_tree_number
##Get one tree from Base_asr.tre and create a file "cafe.tre # head -n3 Base_asr.tre > cafe.tre
tree.cafe<-read.tree("cafe.tre", skip = 2)
tip.order.cafe<-data.frame(node=1:Ntip(tree.cafe), cafe.label = tree.cafe$tip.label)
node.order.cafe<-data.frame(node=1:Nnode(tree.cafe) + Ntip(tree.cafe), cafe.label = tree.cafe$node.label)
nt.order.cafe<-rbind(tip.order.cafe, node.order.cafe)
nt.order.cafe<-separate(nt.order.cafe, cafe.label, into = c("cafe.label"), sep = "_\\d", remove = TRUE)
nt.order.cafe$cafe.label<-gsub("\\*", "", nt.order.cafe$cafe.label)
rm(tree.cafe, tip.order.cafe, node.order.cafe)

##r8s_tree_number
##This is the r8s.ultrametric.tre in CAFE input
tree.r8s<-read.tree("r8s.ultrametric.tre")
myTree<-tree.r8s
tip.order.r8s<-data.frame(node=1:Ntip(tree.r8s), r8s.label = tree.r8s$tip.label)
node.order.r8s<-data.frame(node=1:Nnode(tree.r8s) + Ntip(tree.r8s), r8s.label = tree.r8s$node.label)
nt.order.r8s<-rbind(tip.order.r8s, node.order.r8s)
rm(tree.r8s, tip.order.r8s, node.order.r8s)


##Merge the labels of cafe and orthofinder
nt.order<-merge(nt.order.r8s, nt.order.cafe, by = "node", all.x = TRUE)
rm(nt.order.cafe, nt.order.r8s)

##Duplication number
count.dup<-read.delim("Duplications_per_Species_Tree_Node.tsv", header = TRUE)
names(count.dup)[1]<-"r8s.label"
names(count.dup)[2]<-"Duplications"
names(count.dup)[3]<-"Duplications.5"
nt.order<-merge(nt.order, count.dup, by = "r8s.label", all.x = TRUE)
rm(count.dup)

##Gain and loss number
pvalue=1
pvalue=0.05
ncount=0
ncount=1

##By orthogroups
p.table<-read.delim("Base_family_results.txt", header = TRUE)
p.table<-p.table[p.table$pvalue<pvalue,]
base.change.tab<-read.delim("Base_change.tab", header = FALSE)
names(base.change.tab)<-base.change.tab[1,]
base.change.tab<-base.change.tab[base.change.tab$FamilyID %in% p.table$X.FamilyID, ]
count.gl<-as.data.frame(names(base.change.tab[2:ncol(base.change.tab)]))
names(count.gl)[1]="cafe.label"
count.gl$Increase=NA
count.gl$Decrease=NA
for (i in 1:nrow(count.gl)) {
  count.gl$Increase[i]=nrow(base.change.tab[as.numeric(base.change.tab[,i+1])>(0+ncount),])
  count.gl$Decrease[i]=nrow(base.change.tab[as.numeric(base.change.tab[,i+1])<(0-ncount),])
}

##By each branch
p.table<-read.delim("Base_branch_probabilities.tab", header = FALSE)
names(p.table)<-p.table[1,]
names(p.table)[1]="FamilyID"
p.table<-p.table[p.table$FamilyID != "#FamilyID", ]

base.change.tab<-read.delim("Base_change.tab", header = FALSE)
names(base.change.tab)<-base.change.tab[1,]
base.change.tab<-base.change.tab[base.change.tab$FamilyID != "FamilyID", ]

count.gl<-as.data.frame(names(base.change.tab[2:ncol(base.change.tab)]))
names(count.gl)[1]="cafe.label"
count.gl$Increase=NA
count.gl$Decrease=NA
for (i in 1:nrow(count.gl)) {
  p.table.sub<-p.table[which(as.numeric(p.table[,i+1]) < pvalue),]
  count.gl$Increase[i]=nrow(base.change.tab[as.numeric(base.change.tab[,i+1])>(0+ncount) & base.change.tab$FamilyID %in% p.table.sub$FamilyID,])
  count.gl$Decrease[i]=nrow(base.change.tab[as.numeric(base.change.tab[,i+1])<(0-ncount) & base.change.tab$FamilyID %in% p.table.sub$FamilyID,])
}


#count.gl<-read.delim("Base_clade_results.txt", header = TRUE)
#names(count.gl)[1]<-"cafe.label"
nt.order.new<-merge(nt.order, count.gl, by = "cafe.label", all.x = TRUE)
rm(count.gl)

write.table(nt.order.new, paste0("node_tip_label_count.","ncount",ncount,".pvalue",pvalue,".txt"),
            row.names = FALSE, sep = "\t", 
            quote = FALSE)

##Merge tree
tree<-full_join(myTree, nt.order.new, by = "node")

##No duplication label
{
  p<-ggtree(tree, size = 1, branch.length = "none")+
    geom_tiplab(aes(label=paste("     ",str_replace_all(label,"_"," ")), sep = ""),
                offset = 0.05,
                fontface = "italic",
                size = 6)+
    geom_text2(aes(subset=isTip, label=paste("+",Increase, sep = ""), 
                   hjust = 0.25, vjust = -1),size=4, color = "salmon4")+
    geom_text2(aes(subset=isTip, label=paste("-",Decrease, sep = ""), 
                   hjust = 0.25, vjust = 1.5),size=4, color = "olivedrab4")+
    geom_text2(aes(subset=!isTip & is.na(Duplications)==FALSE, 
                   label=paste("+",Increase, sep = ""), 
                   hjust = 1.1, vjust = -1),
               size=4, color = "salmon4")+
    geom_text2(aes(subset=!isTip & is.na(Duplications)==FALSE, 
                   label=paste("-",Decrease, sep = ""), 
                   hjust = 1.1, vjust = 1.5),
               size=4, color = "olivedrab4")+
    geom_text2(aes(subset=!isTip & is.na(Duplications)==FALSE, 
                   label=r8s.label, 
                   hjust = -0.5, vjust = 0),
               size=4, color = "black")+
    xlim(0, 20)
  
  p
  
  ggsave(paste0("Species_tree_gain_loss_N.","ncount",ncount,".pvalue",pvalue,".png"), 
         width = 10, height = 15, units = "in", dpi = 300)
  ggsave(paste0("Species_tree_gain_loss_N.","ncount",ncount,".pvalue",pvalue,".tiff"), 
         width = 10, height = 15, units = "in", dpi = 300)
  
  f = paste0("Species_tree_gain_loss_N.","ncount",ncount,".pvalue",pvalue,".pptx")
  topptx(p,f, width = 10, height = 15, units = "in")
}

##With duplication
{
  p<-ggtree(tree, size = 1, branch.length = "none")+
    geom_tiplab(aes(label=paste("     ",str_replace_all(label,"_"," ")), sep = ""),
                offset = 0.05,
                fontface = "italic",
                size = 6)+
    geom_label2(aes(subset=isTip, label=Duplications), 
                nudge_x = 0, nudge_y = 0, label.size=0, 
                size=4, color = "black", fill = "lavender")+
    geom_text2(aes(subset=isTip, label=paste("+",Increase, sep = ""), 
                   hjust = 0.2, vjust = -1.2),size=4, color = "salmon4")+
    geom_text2(aes(subset=isTip, label=paste("-",Decrease, sep = ""), 
                   hjust = 0.2, vjust = 2),size=4, color = "olivedrab4")+
    geom_label2(aes(subset=!isTip & is.na(Duplications)==FALSE, label=Duplications), 
                nudge_x = 0, nudge_y = 0, label.size=0, 
                size=4, color = "black", fill = "lavender")+
    geom_text2(aes(subset=!isTip & is.na(Increase)==FALSE, 
                   label=paste("+",Increase, sep = ""), 
                   hjust = 1.1, vjust = -1.2),
               size=4, color = "salmon4")+
    geom_text2(aes(subset=!isTip & is.na(Increase)==FALSE, 
                   label=paste("-",Decrease, sep = ""), 
                   hjust = 1.1, vjust = 2),
               size=4, color = "olivedrab4")+
    geom_text2(aes(subset=!isTip & is.na(Increase)==FALSE, 
                   label=r8s.label, 
                   hjust = -1, vjust = 0),
               size=4, color = "black")+
    xlim(0, 12)
  
  p
  
  ggsave(paste0("Species_tree_gain_loss_dup_N.","ncount",ncount,".pvalue",pvalue,".png"),
         width = 10, height = 20, units = "in", dpi = 300)
  ggsave(paste0("Species_tree_gain_loss_dup_N.","ncount",ncount,".pvalue",pvalue,".tiff"), 
         width = 10, height = 20, units = "in", dpi = 300)
  
  f = paste0(paste0("Species_tree_gain_loss_dup_N.","ncount",ncount,".pvalue",pvalue,".pptx"))
  topptx(p,f, width = 10, height = 20, units = "in")
}

##export plot data
Species_tree_data<-p[["data"]]
write.table(Species_tree_data, file = "Species_tree_gain_loss_dup.txt", row.names = FALSE, sep = "\t")

p<-ggtree(tree, size = 1, branch.length = "none")+
  geom_tiplab(aes(label=paste("     ",str_replace_all(label,"_"," ")), sep = ""),
              offset = 0.05,
              fontface = "italic",
              size = 6, )+
  geom_text2(aes(subset=!isTip & is.na(Duplications)==FALSE, 
                 label=r8s.label, 
                 hjust = -0.2, vjust = 0.5),
             size=6, color = "black")+
  xlim(0, 12)

p

ggsave("Species_tree_N.png", width = 10, height = 15, units = "in", dpi = 300)
ggsave("Species_tree_N.tiff", width = 10, height = 15, units = "in", dpi = 300)

##Ortholog stat
orthostat<-read.delim("Orthogroup.stat.txt")
names(orthostat)[1]="group"
mergestat<-matrix(NA, ncol = 3, nrow = 1, dimnames = list(NA,c("group","species","no")))
mergestat<-as.data.frame(mergestat)
for (i in 2:ncol(orthostat)) {
  substat<-orthostat[,c(1,i)]
  substat$species<-names(substat)[2]
  names(substat)[2]="no"
  mergestat<-rbind(mergestat, substat)
}
mergestat<-subset(mergestat, is.na(no)==FALSE)

speciesorder<-read.delim("Species_tree_gain_loss_dup.txt", header = TRUE)
speciesorder<-subset(speciesorder, node <= ncol(orthostat)-1, select = c("y","label"))
speciesorder<-speciesorder[order(speciesorder$y),2]
#speciesorder<-speciesorder$label
grouporder<-c("Single-copy orthologs",
              "Multi-copy orthologs (>1 species)",
              "Multi-copy orthologs (all species)",
              "Unique paralogs",
              "Unclustered genes")
p<-ggplot(data=mergestat, aes(x=species, y=no, fill=factor(group, levels = grouporder)))+
  geom_bar(stat="identity", width = 0.7)+
  scale_x_discrete(limits = speciesorder, 
                   label = paste("",str_replace_all(speciesorder,"_"," ")))+
  scale_y_continuous(limits = c(0,70000), breaks = seq(0,70000,10000))+
  labs(y="Number of genes", x="")+
  scale_fill_manual(values = c("purple4", "plum2", "burlywood2", "tan4", "darkgray"))+
  theme(axis.line = element_line(linetype = "solid", size = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 0, hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = c(0.75,0.9))

p
ggsave("orthostat.png",
       width = 10, height = 10, units = "in", dpi = 300)
ggsave("orthostat.tiff", 
       width = 10, height = 10, units = "in", dpi = 300)

##Read overlap count table
Orthofinder_dir=""
species.overlaps<-read.table(paste(Orthofinder_dir, 
                                   "Orthogroups_SpeciesOverlaps.tsv", 
                                   sep =""),
                             header = TRUE, sep = "\t")
species.overlaps.1v1<-matrix(data = NA, nrow = 1, ncol = 3,
                             dimnames = list(NA, c("X", "Y","No.overlap.orthogroups")))

for (i in 2:ncol(species.overlaps)) {
  subtable<-species.overlaps[,c(1,i)]
  subtable$Y<-names(subtable)[2]
  names(subtable)[2]="No.overlap.orthogroups"
  species.overlaps.1v1<-rbind(species.overlaps.1v1, subtable)
}
species.overlaps.1v1<-subset(species.overlaps.1v1,
                             is.na(species.overlaps.1v1$X)==FALSE)
species.overlaps.1v1<-unique(species.overlaps.1v1)
rm(subtable)

##Species order
Species_tree_data<-read.table(file = "Species_tree_gain_loss_dup.txt", 
                              header = TRUE, sep = "\t")
nspecies=ncol(species.overlaps)-1
species.order<-subset(Species_tree_data, 
                      node <= nspecies,
                      select = c("y","label"))
species.order<-species.order[order(species.order$y),]$label
species.order
orthogroups<-read.delim("eggnog/Orthogroups.tsv", header = T)
number.orthogroups= nrow(orthogroups)

heatmap.overlap.orthogroups<-ggplot(species.overlaps.1v1,
                                    aes(x=X, y=Y))+
  geom_tile(aes(fill = No.overlap.orthogroups/number.orthogroups))+
  geom_text(aes(label = No.overlap.orthogroups), color = "black", size = 4) +
  scale_fill_gradient(low = "white",high = "mediumpurple4", limits =c(0,0.5))+
  scale_x_discrete(limits = species.order, label=paste("     ",str_replace_all(species.order,"_"," ")))+
  scale_y_discrete(limits = species.order, label=paste("     ",str_replace_all(species.order,"_"," ")))+
  labs(y = "", x = "", fill = "Ratio", title = "")+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, face = "italic", hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, colour = "black", face = "italic"),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.line = element_line(size = 0, colour = NA),
        #panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        panel.background = element_rect(fill = NA))
heatmap.overlap.orthogroups
ggsave("heatmap.overlap.orthogroups.png", width = 16, height = 16, units = "in", dpi = 300)
ggsave("heatmap.overlap.orthogroups.tiff", width = 16, height = 16, units = "in", dpi = 300)
