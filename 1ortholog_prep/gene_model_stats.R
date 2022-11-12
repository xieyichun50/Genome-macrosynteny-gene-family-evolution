##usage
#Rscript gene_model_stats.R -i *.gff3 

##output
#*.gff3.stats.txt

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="standard *.gff3 [default %default]",
              dest="input"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i annotation.gff3 -o out [options]",option_list=option_list)
opt = parse_args(parser)

##manual add files
{
  #opt$input<-"AoH_genomic.gff3"
}

##read in gff3 file
gff3<-read.delim(opt$input, header = FALSE, skip = 1)
names(gff3)<-c("region","source","type",
               "start","stop","score","strand",
               "phase","attributes")

type.list<-unique(gff3$type)
cat(paste0("Type of annotations: ", length(type.list), "\n", "includes: "))
cat(paste0(type.list, collapse = ", "))
cat("\n")

##create the summary table
gff3.summary<-as.data.frame(matrix(data = NA, nrow = 8, ncol = 4))

names(gff3.summary)<-c("Number","N.Average","Length","L.Average")
row.names(gff3.summary)<-c("Gene","Protein coding gene",
                           "mRNA", "Exon", "CDS",
                           "5’-UTR","3’-UTR","tRNA")

##calculation
gff3$length <- gff3$stop-gff3$start+1

####count genes
gff3.summary["Gene","Number"]<-nrow(gff3[gff3$type %in% c("gene","GENE","Gene"),])
gff3.summary["Gene","Length"]<-sum(gff3$length[gff3$type %in% c("gene","GENE","Gene")])
gff3.summary["Gene","L.Average"]<-round(gff3.summary["Gene","Length"]/gff3.summary["Gene","Number"], 1)

####count tRNA
gff3.summary["tRNA","Number"]<-nrow(gff3[gff3$type %in% c("tRNA","TRNA","trna"),])
gff3.summary["tRNA","Length"]<-sum(gff3$length[gff3$type %in% c("tRNA","TRNA","trna")])
gff3.summary["tRNA","L.Average"]<-round(gff3.summary["tRNA","Length"]/gff3.summary["tRNA","Number"], 1)

####count protein coding genes
gff3.summary["Protein coding gene","Number"]<-gff3.summary["Gene","Number"]-gff3.summary["tRNA","Number"]
gff3.summary["Protein coding gene","Length"]<-gff3.summary["Gene","Length"]-gff3.summary["tRNA","Length"]
gff3.summary["Protein coding gene","L.Average"]<-round(gff3.summary["Protein coding gene","Length"]/
                                                         gff3.summary["Protein coding gene","Number"], 1)

####count mRNA of protein coding genes
gff3.summary["mRNA","Number"]<-nrow(gff3[gff3$type %in% c("mRNA","MRNA","mrna"),])
gff3.summary["mRNA","N.Average"]<-round(gff3.summary["mRNA","Number"]/gff3.summary["Protein coding gene","Number"],2)
gff3.summary["mRNA","Length"]<-sum(gff3$length[gff3$type %in% c("mRNA","MRNA","mrna")])
gff3.summary["mRNA","L.Average"]<-round(gff3.summary["mRNA","Length"]/gff3.summary["mRNA","Number"], 1)

####count exon of protein coding genes
gff3.summary["Exon","Number"]<-nrow(gff3[gff3$type %in% c("Exon","EXON","exon"),])-gff3.summary["tRNA","Number"]
gff3.summary["Exon","N.Average"]<-round(gff3.summary["Exon","Number"]/gff3.summary["mRNA","Number"],2)
gff3.summary["Exon","Length"]<-sum(gff3$length[gff3$type %in% c("Exon","EXON","exon")])-gff3.summary["tRNA","Length"]
gff3.summary["Exon","L.Average"]<-round(gff3.summary["Exon","Length"]/gff3.summary["Exon","Number"], 1)

####count cds of protein coding genes
gff3.summary["CDS","Number"]<-nrow(gff3[gff3$type %in% c("cds","CDS"),])
gff3.summary["CDS","N.Average"]<-round(gff3.summary["CDS","Number"]/gff3.summary["mRNA","Number"],2)
gff3.summary["CDS","Length"]<-sum(gff3$length[gff3$type %in% c("cds","CDS")])
gff3.summary["CDS","L.Average"]<-round(gff3.summary["CDS","Length"]/gff3.summary["CDS","Number"], 1)

####count 5'-UTR of protein coding genes
gff3.summary["5’-UTR","Number"]<-nrow(gff3[gff3$type %in% 
                                             c("five_prime_UTR","Five_prime_UTR","5-UTR","5-utr","5’-UTR","5’-utr"),])
gff3.summary["5’-UTR","N.Average"]<-round(gff3.summary["5’-UTR","Number"]/gff3.summary["mRNA","Number"],2)
gff3.summary["5’-UTR","Length"]<-sum(gff3$length[gff3$type %in% 
                                                   c("five_prime_UTR","Five_prime_UTR","5-UTR","5-utr","5’-UTR","5’-utr")])
gff3.summary["5’-UTR","L.Average"]<-round(gff3.summary["5’-UTR","Length"]/gff3.summary["5’-UTR","Number"], 1)

####count 3'-UTR of protein coding genes
gff3.summary["3’-UTR","Number"]<-nrow(gff3[gff3$type %in% 
                                             c("three_prime_UTR","Three_prime_UTR",
                                               "3-UTR","3-utr","3’-UTR","3’-utr"),])
gff3.summary["3’-UTR","N.Average"]<-round(gff3.summary["3’-UTR","Number"]/gff3.summary["mRNA","Number"],2)
gff3.summary["3’-UTR","Length"]<-sum(gff3$length[gff3$type %in% 
                                                   c("three_prime_UTR","Three_prime_UTR",
                                                     "3-UTR","3-utr","3’-UTR","3’-utr")])
gff3.summary["3’-UTR","L.Average"]<-round(gff3.summary["3’-UTR","Length"]/gff3.summary["3’-UTR","Number"], 1)

cat(paste0("File will be written to ",opt$input,".stats.txt", "\n"))

write.table(gff3.summary, paste0(opt$input,".stats.txt"),
            quote = F, sep = "\t")
cat(paste0("Done! \n"))
