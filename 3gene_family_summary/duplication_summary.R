##Summarise the duplication count table
dup<-read.delim("F:/Snail/longest_protO/Orthofinder/Results_May27/Gene_Duplication_Events/Duplications.tsv")
dup.summary<-xtabs(~ Orthogroup+Species.Tree.Node, data = dup)
dup.summary<-as.data.frame(dup.summary)
dup.table<-matrix(NA, 
                  nrow = length(unique(dup.summary$Orthogroup)), 
                  ncol = 1,
                  dimnames = list(unique(dup.summary$Orthogroup), "Orthogroup"))
dup.table<-as.data.frame(dup.table)
dup.table$Orthogroup<-unique(dup.summary$Orthogroup)
nodelist<-as.character(unique(dup.summary$Species.Tree.Node))
nodelist<-as.data.frame(nodelist)
for (i in 1:nrow(nodelist)) {
  dupsub<-subset(dup.summary, 
                 Species.Tree.Node == nodelist$nodelist[i], select = c("Orthogroup", "Freq"))
  names(dupsub)[2]<-nodelist$nodelist[i]
  dup.table<-merge(dup.table, dupsub, by = "Orthogroup", all.x = TRUE)
}
write.table(dup.table,
            "duplication_count.txt",
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)
