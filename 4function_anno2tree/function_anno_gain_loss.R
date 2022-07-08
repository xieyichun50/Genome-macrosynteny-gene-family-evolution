library(dplyr)
library(tidyr)
library(eoffice)

Genesorthogrouppair.1v1<-read.delim("all.Genesorthogrouppair.1v1", header = TRUE, sep = "\t")
GenesGOpair.1v1<-read.delim("all.GenesGOpair.1v1", header = TRUE, sep = "\t")
GenesKOGpair.1v1<-read.delim("all.GenesKOGpair.1v1", header = TRUE, sep = "\t")
GenesKEGGpair.1v1<-read.delim("all.GenesKEGGpair.1v1", header = TRUE, sep = "\t")
#Read change table
#change.count<-read.delim("/home/yichun/huilab/Eurema_hecabe/longest_prot/OrthoFinder/Results_May05/CAFE/r8s_lambda3/Gamma_change.tab", header = TRUE)
change.count<-read.delim("Gamma_change.tab", header = TRUE)

names(change.count)[1]="orthogroups"
change.count.KOG<-merge(change.count, OG.KOGpair, by = "orthogroups", all.x = TRUE)
change.count.KOG<-subset(change.count.KOG, is.na(KOG)==FALSE)

kog2name<-read.delim("/home/yichun/tools/Annotation/kog2name.txt",
                     sep = "\t", colClasses = "character")
KOG.change.count<-rbind(kog2name, kog2name)
KOG.change.count$KOG[1:25]<-paste(KOG.change.count$kogClass[1:25], ".gain", sep = "")
KOG.change.count$KOG[26:50]<-paste(KOG.change.count$kogClass[26:50], ".loss", sep = "")

for (i in 2:(ncol(change.count.KOG)-1)) {
  change.count.KOG.sub<-change.count.KOG[,c(1, i, ncol(change.count.KOG))]
  a=names(change.count.KOG.sub)[2]
  names(change.count.KOG.sub)[2]="count"
  change.count.KOG.sub$count<-as.numeric(change.count.KOG.sub$count)

  change.count.KOG.subgain<-change.count.KOG.sub[which(change.count.KOG.sub$count>0),]
  change.count.KOG.subloss<-change.count.KOG.sub[which(change.count.KOG.sub$count<0),]
  KOG.change.count.subgain<-change.count.KOG.subgain %>% group_by(KOG) %>% summarise(sum = sum(count))
  KOG.change.count.subloss<-change.count.KOG.subloss %>% group_by(KOG) %>% summarise(sum = sum(count))
  KOG.change.count.subgain$KOG<-paste(KOG.change.count.subgain$KOG, ".gain", sep = "")
  names(KOG.change.count.subgain)[2]=a
  KOG.change.count.subloss$KOG<-paste(KOG.change.count.subloss$KOG, ".loss", sep = "")
  names(KOG.change.count.subloss)[2]=a
  KOG.change.count.sub<-rbind(KOG.change.count.subgain, KOG.change.count.subloss)
  KOG.change.count<-merge(KOG.change.count, KOG.change.count.sub, by = "KOG", all.x = TRUE)
}
rm(change.count.KOG.sub, change.count.KOG.subgain, change.count.KOG.subloss,
   KOG.change.count.subgain, KOG.change.count.subloss, KOG.change.count.sub)

write.table(KOG.change.count, file = "KOG.change.count", row.names = FALSE, sep = "\t")
