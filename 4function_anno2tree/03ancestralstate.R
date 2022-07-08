##usage
#Rscript /jelly_data/yichun/scripts/5formatenrich/03ancestralstate.R -c Base_count.tab -d Base_change.tab -g Orthogroups.GeneCount.tsv -n Species_tree_gain_loss_dup.txt
#Output the following
##pseudo.species.txt --pseudogenome of each node/tip of tree
##gain.species.txt --gained genes of each node/tip of tree
##loss.species.txt --lost genes of each node/tip of tree
##genome.species.txt --genomes of species
##specific.species.txt --specific orthologues of species

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-c","--count"), type="character", default=NULL,
              help="base_count.tab from cafe5' [default %default]",
              dest="count.table"),
  make_option(c("-d","--delta"), type="character", default=NULL,
              help="base_change.tab from cafe5' [default %default]",
              dest="change.table"),
  make_option(c("-g","--genome"), type="character", default=NULL,
              help="Orthogroups.GeneCount.tsv from orthofinder' [default %default]",
              dest="genome.table"),
  make_option(c("-n","--nodelabel"), type="character", default=NULL,
              help="Species_tree_gain_loss_dup.txt from gain loss tree' [default %default]",
              dest="speciesorder"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##Manually add files
{
  #setwd("D:/centipede/cafe")
  #opt$count.table="Base_count.tab"
  #opt$change.table="Base_change.tab"
  #opt$genome.table="Orthogroups.GeneCount.tsv"
  #opt$speciesorder="Species_tree_gain_loss_dup.txt"
}

#read in species label
speciesorder<-read.delim(opt$speciesorder, header = TRUE, stringsAsFactors = FALSE)

#read in count table and format node labels
{
  count.table<-read.delim(opt$count.table, header = FALSE, stringsAsFactors = FALSE)
  row.names(count.table)<-count.table$V1
  names(count.table)<-count.table[1,]
  count.table<-count.table[-1,-1]
  for (i in 1:nrow(speciesorder)) {
    names(count.table)[names(count.table)==speciesorder$cafe.label[i]]=speciesorder$label[i]
  }
  
  count.table<-data.frame(lapply(count.table,as.numeric), row.names = row.names(count.table))
}

#read in change table and format node labels
{
  change.table<-read.delim(opt$change.table, header = FALSE, stringsAsFactors = FALSE)
  row.names(change.table)<-change.table$V1
  names(change.table)<-change.table[1,]
  change.table<-change.table[-1,-1]
  for (i in 1:nrow(speciesorder)) {
    names(change.table)[names(change.table)==speciesorder$cafe.label[i]]=speciesorder$label[i]
  }
  
  change.table<-data.frame(lapply(change.table,as.numeric), row.names = row.names(change.table))
}

#read in genome count table and format node labels
{
  genome.table<-read.delim(opt$genome.table, header = TRUE, stringsAsFactors = FALSE)
  row.names(genome.table)<-genome.table$Orthogroup
  genome.table<-genome.table[,-1]
  genome.table<-data.frame(lapply(genome.table,as.numeric), row.names = row.names(genome.table))
  
}

##Test
##Are loss number < parent node gene number, should be all 0
{
  for (j in 1:ncol(count.table)){
    cat(paste0(names(count.table)[j], " = "))
    parentnode<-subset(speciesorder, label == names(count.table)[j])
    parentspecies<-speciesorder$label[which(speciesorder$node==parentnode$parent)]
    testcount<-count.table[,which(names(count.table)==parentspecies)]+change.table[,which(names(change.table)==names(count.table)[j])]
    cat(length(which(testcount<0)==TRUE,"\n"))
  }
}
##Are gain number < current node gene number, should be all 0
{
  for (j in 1:ncol(count.table)){
    cat(paste0(names(count.table)[j], " = "))
    testcount<-count.table[,j]-change.table[,which(names(change.table)==names(count.table)[j])]
    cat(length(which(testcount<0)==TRUE),"\n")
  }
}

#Generate ancestral state
##Create blank matrix
{
  pseudogenome<-matrix(NA, nrow = 1, ncol = 2,
                       dimnames = list(NA,c("Genes","Orthogroup")))
  pseudogenome<-as.data.frame(pseudogenome)
  #pseudogenome<-pseudogenome[which(is.na(pseudogenome$Genes)==FALSE),]
  
  pseudo.all<-pseudogenome
  gain.all<-pseudogenome
  loss.all<-pseudogenome
}

##Main pool ##No need
{
  for (i in 1:nrow(count.table)) {
    cat(paste0("Generating pseudogenome of all", "\n"))
    rowmax<-max(count.table[i,1:ncol(count.table)])
    OG.sub<-matrix(paste(row.names(count.table)[i],".",c(1:rowmax), sep = ""))
    OG.sub<-as.data.frame(OG.sub)
    names(OG.sub)[1]<-"Genes"
    OG.sub$Orthogroup<-row.names(count.table)[i]
    pseudo.all<-rbind(pseudo.all, OG.sub)
  }
  pseudo.all<-pseudo.all[which(is.na(pseudo.all$Genes)==FALSE),]
  write.table(pseudo.all, 
              file = "pseudo.all.txt",
              row.names = FALSE, sep = "\t", quote = FALSE)
}

##Genome pools
{
  for (j in 1:ncol(count.table)){
    cat(paste0("Generating pseudogenome of ", names(count.table)[j], "\n"))
    pseudo.species.genes<-pseudogenome
    count.table.sub<-data.frame(count.table[,j], row.names = row.names(count.table))
    speciesname=names(count.table)[j]
    names(count.table.sub)[1]="Orthogroup"
    count.table.sub<-subset(count.table.sub, Orthogroup > 0)
    for (i in 1:nrow(count.table.sub)) {
      n<-count.table.sub[i,1]
      OG.sub<-matrix(paste(row.names(count.table.sub)[i],".",c(1:n), sep = ""))
      OG.sub<-as.data.frame(OG.sub)
      names(OG.sub)[1]<-"Genes"
      OG.sub$Orthogroup<-row.names(count.table.sub)[i]
      pseudo.species.genes<-rbind(pseudo.species.genes, OG.sub)
    }
    pseudo.species.genes<-pseudo.species.genes[which(is.na(pseudo.species.genes$Genes)==FALSE),]
    write.table(pseudo.species.genes, 
                file = paste("pseudo.",speciesname,".txt", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
    rm(pseudo.species.genes, count.table.sub, speciesname, n, OG.sub)
  }
}

##Gain loss pool
{
  for (j in 1:ncol(change.table)){
    cat(paste0("Generating gain and loss lists of ", names(count.table)[j], "\n"))
    gain.species.genes<-pseudogenome
    loss.species.genes<-pseudogenome
    change.table.sub<-data.frame(change.table[,j], row.names = row.names(change.table))
    speciesname=names(change.table)[j]
    names(change.table.sub)[1]="Orthogroup"
    change.table.sub<-subset(change.table.sub, Orthogroup != 0)
    for (i in 1:nrow(change.table.sub)) {
      n<-change.table.sub[i,1]
      OG.sub<-matrix(paste(row.names(change.table.sub)[i],".",c(1:abs(n)), sep = ""))
      OG.sub<-as.data.frame(OG.sub)
      names(OG.sub)[1]<-"Genes"
      OG.sub$Orthogroup<-row.names(change.table.sub)[i]
      if (n>0) {
        gain.species.genes<-rbind(gain.species.genes, OG.sub)
      } else if (n<0) {
        loss.species.genes<-rbind(loss.species.genes, OG.sub)
      } else {
      }
    }
    gain.species.genes<-gain.species.genes[which(is.na(gain.species.genes$Genes)==FALSE),]
    loss.species.genes<-loss.species.genes[which(is.na(loss.species.genes$Genes)==FALSE),]
    write.table(gain.species.genes, 
                file = paste("gain.",speciesname,".txt", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(loss.species.genes, 
                file = paste("loss.",speciesname,".txt", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
    rm(gain.species.genes, loss.species.genes, change.table.sub, speciesname, n, OG.sub)
  }
}

##Species pool
{
  for (j in 1:(ncol(genome.table)-1)){
    cat(paste0("Generating species genelist and species-specific genelist of ", names(count.table)[j], "\n"))
    genome.species.genes<-pseudogenome
    specific.species.genes<-pseudogenome
    genome.table.sub<-data.frame(genome.table[,c(j, ncol(genome.table))], row.names = row.names(genome.table))
    speciesname=names(genome.table)[j]
    names(genome.table.sub)[1]="Orthogroup"
    genome.table.sub<-subset(genome.table.sub, Orthogroup > 0)
    for (i in 1:nrow(genome.table.sub)) {
      n<-genome.table.sub[i,1]
      OG.sub<-matrix(paste(row.names(genome.table.sub)[i],".",c(1:n), sep = ""))
      OG.sub<-as.data.frame(OG.sub)
      names(OG.sub)[1]<-"Genes"
      OG.sub$Orthogroup<-row.names(genome.table.sub)[i]
      genome.species.genes<-rbind(genome.species.genes, OG.sub)
      sc<-genome.table.sub$Total[i]-genome.table.sub[i,1]
      if (sc == 0) {
        specific.species.genes<-rbind(specific.species.genes, OG.sub)
      } else {
        
      }
    }
    genome.species.genes<-genome.species.genes[which(is.na(genome.species.genes$Genes)==FALSE),]
    specific.species.genes<-specific.species.genes[which(is.na(specific.species.genes$Genes)==FALSE),]
    write.table(genome.species.genes, 
                file = paste("genome.",speciesname,".txt", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(specific.species.genes, 
                file = paste("specific.",speciesname,".txt", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
    rm(genome.species.genes, specific.species.genes, genome.table.sub, speciesname, n, OG.sub)
  }
}

