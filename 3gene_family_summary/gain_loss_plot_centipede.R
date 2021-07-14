library(ape)
library(ggtree)
library(flextable)
setwd("D:/2gain_loss/centipede")

myTree<-read.tree(file = "cafe.number.node.tree")
gl<-read.delim("Gamma_clade_results.txt")
gl.node<-gl
names(gl.node)[1]="cafe.label"
gl.tip<-gl
names(gl.tip)[1]="cafe.label"
node.order<-data.frame(node=1:Nnode(myTree) + Ntip(myTree), cafe.label = myTree$node.label)
node.order<-merge(node.order, gl.node, by = "cafe.label", all.x = TRUE)
tip.order<-data.frame(node=1:Ntip(myTree), cafe.label = myTree$tip.label)
tip.order<-merge(tip.order, gl.tip, by = "cafe.label", all.x = TRUE)
nt.order<-rbind(node.order, tip.order)
nt.order<-nt.order[,2:4]
tree<-full_join(myTree, nt.order, by = "node")


p<-ggtree(tree)+
  geom_tiplab(size=3, fontface = "italic")+
  geom_text2(aes(subset=isTip, label=paste("+",Increase,"/-",Decrease, sep = ""), hjust = 0, vjust = 2),size=3, color = "blue")+
  geom_text2(aes(subset=(!isTip & node !=10 & node !=11 & node != 16), label=paste("+",Increase,"/-",Decrease, sep = ""), hjust = 1, vjust = -1),size=3, color = "blue")+
  geom_text2(aes(subset=node==11, label=paste("+",Increase,"/-",Decrease, sep = ""), hjust = 1, vjust = 2),size=3, color = "blue")+
  geom_text2(aes(subset=node==16, label=paste("+",Increase,"/-",Decrease, sep = ""), hjust = 1, vjust = 2),size=3, color = "blue")+
  geom_text2(aes(subset=node==10, label=paste("+",Increase,"/-",Decrease, sep = ""), hjust = 1, vjust = 1),size=3, color = "NA")+
  
  xlim(-0.2, 1.2)
p
ggsave("Species_tree_gain_loss.png", width = 8, height = 6, units = "in", dpi = 300)

