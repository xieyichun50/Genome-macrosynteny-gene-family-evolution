#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="all-vs-all.blast result to two species' [default %default]",
              dest="input_filename"),
  make_option(c("-o","--output"), type="character", default="out",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]"),
  make_option(c("-A", "--xspecies-fai"), type="character", default=NULL,
              help="fasta index of reference (x) genome [default %default]",
              dest="fai_ref"),
  make_option(c("-B", "--yspecies-fai"), type="character", default=NULL,
              help="fasta index of query (y) genome [default %default]",
              dest="fai_query"),
  make_option(c("-r", "--xspecies-gff3"), type="character", default=NULL,
              help="mRNA gff3 of reference (x) genome [default %default]",
              dest="gff3_ref"),
  make_option(c("-q", "--yspecies-gff3"), type="character", default=NULL,
              help="mRNA gff3 of query (y) genome [default %default]",
              dest="gff3_query"),
  make_option(c("-R", "--xspecies-locustag"), type="character", default=NULL,
              help="locustag of reference (x) genome [default %default]",
              dest="locustag_ref"),
  make_option(c("-Q", "--yspecies-locustag"), type="character", default=NULL,
              help="locustag of query (y) genome [default %default]",
              dest="locustag_query"),
  make_option(c("-a", "--xspecies-name"), type="character", default="x species",
              help="Reference name show at xlab [default %default]",
              dest="name_ref"),
  make_option(c("-b", "--yspecies-name"), type="character", default="y species",
              help="query name show at ylab [default %default]",
              dest="name_qry"),
  make_option(c("-m", "--min-orthologs-x"), type="numeric", default=20,
              help="minimum orthologs number in reference scaffolds [default %default]",
              dest="min_ortho_x"),
  make_option(c("-n", "--min-orthologs-y"), type="numeric", default=5,
              help="minimum orthologs number in query scaffolds [default %default]",
              dest="min_ortho_y"),
  make_option(c("-p","--plot-size"), type="numeric", default=15,
              help="plot size X by X inches [default %default]",
              dest="plot_size"),
  make_option(c("-x", "--show-horizontal-lines"), action="store_true", default=TRUE,
              help="turn on horizontal lines on plot for separating scaffolds  [default %default]",
              dest="h_lines"),
  make_option(c("-y", "--show-verticle-lines"), action="store_true", default=TRUE,
              help="turn on verticle lines on plot for separating scaffolds  [default %default]",
              dest="v_lines"),
  make_option(c("-z", "--break-size"), type="numeric", default=50,
              help="second y break by z Mb [default %default]",
              dest="break_size"),
  make_option(c("-w", "--window-size"), type="numeric", default=50,
              help="sliding window size for significance test [default %default]",
              dest="window")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

setwd("/home/yichun/kairuku/huilab/synteny202103/longest_prot/OrthoFinder/Results_Mar19/Orthologues/formatted")
opt$input_filename="ortholog.Lithobius_niger__v__Glomeris_maerens.tsv.txt"
opt$fai_ref="/home/yichun/huilab/centipede/lni.fa.fai"
opt$fai_query="/home/yichun/huilab/centipede/gma.fa.fai"
opt$gff3_ref="/home/yichun/huilab/centipede/lni.mRNA.gff3"
opt$gff3_query="/home/yichun/huilab/centipede/gma.mRNA.gff3"
opt$name_ref="Lithobius niger"
opt$name_qry="Glomeris maerens"
opt$locustag_ref="Lni"
opt$locustag_query="Gma"
opt$plot_size=15
opt$break_size=50
opt$min_ortho_x=5
#opt$window=50
#opt$output_filename="test"

##print parameters to terminal
cat(paste0("PARAMETERS:\ninput (-i): ", opt$input_filename,"\n"))
cat(paste0("output (-o): ", opt$output_filename,"\n"))
cat(paste0("xspecies-fai (-A): ", opt$fai_ref,"\n"))
cat(paste0("yspecies-fai (-B): ", opt$fai_query,"\n"))
cat(paste0("xspecies-gff3 (-r): ", opt$gff3_ref,"\n"))
cat(paste0("yspecies-gff3 (-q): ", opt$gff3_query,"\n"))
cat(paste0("xspecies-locustag (-R): ", opt$locustag_ref,"\n"))
cat(paste0("yspecies-locustag (-Q): ", opt$locustag_query,"\n"))
cat(paste0("xspecies-name (-a): ", opt$name_ref,"\n"))
cat(paste0("yspecies-name (-b): ", opt$name_qry,"\n"))
cat(paste0("minimum orthologs number in reference scaffolds (-m): ", opt$min_ortho_x,"\n"))
cat(paste0("minimum orthologs number in query scaffolds (-n): ", opt$min_ortho_y,"\n"))
cat(paste0("sliding window size for significance test (-w): ", opt$window,"\n"))
cat(paste0("plot size (-p): ", opt$plot_size,"\n"))
cat(paste0("show-horizontal-lines (-x): ", opt$h_lines,"\n"))
cat(paste0("show-verticle-lines (-y): ", opt$v_lines,"\n"))

##Ortholog file formatting
orthologs<-read.table(opt$input_filename, stringsAsFactors = F)
#names(orthologs)[3]="score"
#orthologs<-separate(orthologs, "V1", c("x.sp", "x.gene"),
#                    sep = "\\|", remove = TRUE, convert = TRUE)
#orthologs<-separate(orthologs, "V2", c("y.sp", "y.gene"),
#                    sep = "\\|", remove = TRUE, convert = TRUE)

#orthologs<-subset(orthologs, select = c("x.gene","y.gene","score"))

names(orthologs)[1]="x.gene"
names(orthologs)[2]="y.gene"
####in case of trans and reverse, copy and merge
orthologs.reverse<-subset(orthologs, select = c("y.gene","x.gene"))
names(orthologs.reverse)[1]="x.gene"
names(orthologs.reverse)[2]="y.gene"

orthologs.merge<-rbind(orthologs,orthologs.reverse)
rm(orthologs, orthologs.reverse)

####get only orthologs of targeted species
orthologs.subset<-orthologs.merge[grep(pattern=opt$locustag_ref,orthologs.merge$x.gene),]
orthologs.subset<-orthologs.subset[grep(pattern=opt$locustag_query,orthologs.subset$y.gene),]
rm(orthologs.merge)
##Ortholog file formatted

##gff3 files formatting 
####reference genome gff3
gff3_ref<-read.table(opt$gff3_ref, stringsAsFactors = F)
gff3_ref<-gff3_ref[,c(1,4,5,9)]
gff3_ref<-separate(gff3_ref, "V9", c("x.gene","other"),
                   sep = ";", remove = TRUE, convert = TRUE)
gff3_ref<-separate(gff3_ref, "x.gene", c("p","x.gene"),
                   sep = "=", remove = TRUE, convert = TRUE)
names(gff3_ref)[1]="x.chr"
names(gff3_ref)[2]="x.start"
names(gff3_ref)[3]="x.end"

####for (x,y) position
gff3_ref<-gff3_ref[,c("x.chr", "x.start", "x.end", "x.gene")]

####query genome gff3
gff3_query<-read.table(opt$gff3_query, stringsAsFactors = F)
gff3_query<-gff3_query[,c(1,4,5,9)]
gff3_query<-separate(gff3_query, "V9", c("y.gene","other"),
                     sep = ";", remove = TRUE, convert = TRUE)
gff3_query<-separate(gff3_query, "y.gene", c("p","y.gene"),
                     sep = "=", remove = TRUE, convert = TRUE)
names(gff3_query)[1]="y.chr"
names(gff3_query)[2]="y.start"
names(gff3_query)[3]="y.end"

####for (x,y) position
gff3_query<-gff3_query[,c("y.chr", "y.start", "y.end", "y.gene")]

##gff3 files formatted

##Merge orthologs with position and filter out orthologs on scaffolds with a few orthologs
orthologs.subset<-merge(orthologs.subset, gff3_ref, by = "x.gene", all.x = TRUE)
orthologs.subset<-merge(orthologs.subset, gff3_query, by = "y.gene", all.x = TRUE)

x.chrID<-as.list(xtabs(~x.chr, data = orthologs.subset))
x.chrID<-x.chrID[x.chrID>=opt$min_ortho_x] ##keep only x scaffolds with opt$min_ortho_x orthologs
orthologs.filtered = orthologs.subset[which(orthologs.subset$x.chr %in% names(x.chrID)),]

y.chrID<-as.list(xtabs(~y.chr, data = orthologs.filtered))
y.chrID<-y.chrID[y.chrID>=opt$min_ortho_y] ##keep only y scaffolds with min_ortho_y orthologs
orthologs.filtered = orthologs.filtered[which(orthologs.filtered$y.chr %in% names(y.chrID)),]
orthologs.filtered = orthologs.filtered[order(orthologs.filtered$x.start, orthologs.filtered$x.chr, 
                                              decreasing = FALSE),]
##ortholog filter finished

##Import chr/scaffold length and renumber
fai_ref<-read.table(opt$fai_ref, stringsAsFactors = F)
fai_ref<-fai_ref[,c(1,2)]
names(fai_ref)[1]="x.chr"
names(fai_ref)[2]="x.length"
fai_ref<-fai_ref[which(fai_ref$x.chr %in% names(x.chrID)),]
fai_ref<-fai_ref[order(fai_ref$x.length, decreasing = TRUE),]
fai_ref$x.rank<-1:nrow(fai_ref)
row.names(fai_ref)<-1:nrow(fai_ref)
fai_ref$x.chrstart=NA
fai_ref$x.chrstart[1]<-0
fai_ref$x.chrend=NA
fai_ref$x.chrend[1]=fai_ref$x.length[1]
for (i in 2:nrow(fai_ref)) {
  fai_ref$x.chrstart[i]<-fai_ref$x.length[i-1]+fai_ref$x.chrstart[i-1]
  fai_ref$x.chrend[i]<-fai_ref$x.length[i]+fai_ref$x.chrend[i-1]
}
max(fai_ref$x.chrend)

fai_query<-read.table(opt$fai_query, stringsAsFactors = F)
fai_query<-fai_query[,c(1,2)]
names(fai_query)[1]="y.chr"
names(fai_query)[2]="y.length"
fai_query<-fai_query[which(fai_query$y.chr %in% names(y.chrID)),]

####order y axis
chr.freq<-as.data.frame(xtabs(~x.chr+y.chr, data = orthologs.filtered))
chr.freq<-chr.freq[order(chr.freq$y.chr,chr.freq$Freq, decreasing = TRUE),]
chr.freq<-chr.freq[which(chr.freq$x.chr %in% fai_ref$x.chr),]
row.names(chr.freq)<-1:nrow(chr.freq)

fai_query$x.chr<-NA
for (i in 1:nrow(fai_query)) {
  a<-subset(chr.freq, chr.freq$y.chr==fai_query$y.chr[i])
  a<-a[order(a$Freq, decreasing = TRUE),]
  row.names(a)<-1:nrow(a)
  fai_query$x.chr[i]<-as.character(a[1,1])
}
fai_query<-merge(fai_query, fai_ref[,c("x.chr","x.rank")], by = "x.chr", all.x = TRUE)
fai_query$y.rank<-NA
for (i in 1:nrow(fai_query)) {
  if (is.na(fai_query$y.rank[i])==TRUE) {
    b<-subset(chr.freq, chr.freq$x.chr==fai_query$x.chr[i] & chr.freq$Freq > 0)
    b<-b[order(b$Freq, decreasing = TRUE),]
    b$y.rate<-1:nrow(b)
    fai_query<-merge(fai_query, b[,c("x.chr", "y.chr", "y.rate")], by = c("x.chr","y.chr"), all.x = TRUE)
    if (is.na(fai_query$y.rank[i])==TRUE) {
      fai_query$y.rank[i]<-fai_query$y.rate[i]
      fai_query<-fai_query[,c(1:5)]
    } else {}
    
  } else{}
  
}

fai_query<-fai_query[order(fai_query$x.rank, fai_query$y.rank, decreasing = FALSE),]
#fai_query<-fai_query[order(fai_query$y.length, decreasing = TRUE),]
row.names(fai_query)<-1:nrow(fai_query)
fai_query$y.chrstart=NA
fai_query$y.chrstart[1]<-0
fai_query$y.chrend=NA
fai_query$y.chrend[1]=fai_query$y.length[1]
for (i in 2:nrow(fai_query)) {
  fai_query$y.chrstart[i]<-fai_query$y.length[i-1]+fai_query$y.chrstart[i-1]
  fai_query$y.chrend[i]<-fai_query$y.length[i]+fai_query$y.chrend[i-1]
}
max(fai_query$y.chrend)
##Chr/scaffold relengthed

####Significance test
chr.windows.p<-matrix(data = NA, ncol = 7, 
                      dimnames = list(c(), c("x.chr","y.chr","a","b","c","d","pvalue")))
chr.windows.p<-as.data.frame(chr.windows.p)
p.raw<-chr.windows.p
n.ortholog<-nrow(orthologs.filtered)
for (i in 1:nrow(fai_ref)) {
  for (j in 1:nrow(fai_query)) {
    pool.ortholog<-subset(orthologs.filtered, 
                          x.chr == fai_ref$x.chr[i] & y.chr == fai_query$y.chr[j])
    ####not count chr pairs with 0 ortholog
    if (nrow(pool.ortholog > 0)) {
      row.names(pool.ortholog)<-1:nrow(pool.ortholog)
      pool.x.genes<-subset(orthologs.filtered, x.chr == fai_ref$x.chr[i])
      pool.x.genes<-pool.x.genes[order(pool.x.genes$x.start),]
      row.names(pool.x.genes)<-1:nrow(pool.x.genes)
      n.pool.x.genes<-nrow(pool.x.genes)
      
      pool.y.genes<-subset(orthologs.filtered, y.chr == fai_query$y.chr[j])
      pool.y.genes<-pool.y.genes[order(pool.y.genes$y.start),]
      row.names(pool.y.genes)<-1:nrow(pool.y.genes)
      n.pool.y.genes<-nrow(pool.y.genes)
      ####chr.x and chr.y all not be split
      pool.mn.genes<-pool.ortholog
      a<-nrow(pool.mn.genes)
      b<-n.pool.x.genes-a
      c<-n.pool.y.genes-a
      d<-n.ortholog-a-b-c
      p<-p.raw
      p[1,1]=fai_ref$x.chr[i]
      p[1,2]=fai_query$y.chr[j]
      p[1,3]=a
      p[1,4]=b
      p[1,5]=c
      p[1,6]=d
      t<-fisher.test(matrix(c(a,b,c,d), ncol = 2, nrow = 2), alternative = "greater")
      p[1,7]=t$p.value
      chr.windows.p<-rbind(chr.windows.p,p)
    } else {}
  }
}

####P value computed
chr.windows.p$pvalue<-as.numeric(chr.windows.p$pvalue)
####adjust P value with BH method
chr.windows.p$pvalue.adjust<-p.adjust(chr.windows.p$pvalue, method = "BH", n = nrow(chr.windows.p))
####significant genes
orthologs.filtered<-merge(orthologs.filtered, 
                          chr.windows.p[,c("x.chr","y.chr","pvalue.adjust")],
                          by = c("x.chr","y.chr"),
                          all.x = TRUE)

##Assign new position to orthologs
####Merge x-reference
orthologs.filtered<-merge(orthologs.filtered, fai_ref[,c("x.chr","x.chrstart")], by = "x.chr", all.x = TRUE)
orthologs.filtered$x.pos.new<-(orthologs.filtered$x.start+orthologs.filtered$x.end)/2+orthologs.filtered$x.chrstart

####Merge y-query
orthologs.filtered<-merge(orthologs.filtered, fai_query[,c("y.chr","y.chrstart")], by = "y.chr", all.x = TRUE)
orthologs.filtered$y.pos.new<-(orthologs.filtered$y.start+orthologs.filtered$y.end)/2+orthologs.filtered$y.chrstart
##New (x,y) generated

##Plotting
p<-ggplot()+
  geom_point(data = subset(orthologs.filtered, pvalue.adjust < 0.05), 
             aes(x=x.pos.new/10^6, 
                 y=y.pos.new/10^6, 
                 color = factor(x.chr)), size = 2)+
  geom_point(data = subset(orthologs.filtered, pvalue.adjust >= 0.05), 
             aes(x=x.pos.new/10^6, 
                 y=y.pos.new/10^6), size = 2, color = "Gray50")+
  scale_x_continuous(breaks = fai_ref$x.chrend/10^6,
                     labels = fai_ref$x.chr,
                     sec.axis = sec_axis(~., name = "Mb",
                                         breaks = seq(0, 
                                                      ceiling(max(fai_ref$x.chrend)/10^6), 
                                                      ceiling(max(fai_ref$x.chrend)/(10^6*8*50))*50)))+
  scale_y_continuous(breaks = fai_query$y.chrend/10^6,
                     labels = fai_query$y.chr,
                     sec.axis = sec_axis(~., name = "Mb",
                                         breaks = seq(0, 
                                                      ceiling(max(fai_query$y.chrend)/10^6), 
                                                      ceiling(max(fai_query$y.chrend)/(10^6*8*opt$break_size))*opt$break_size)))+
  labs(title = paste0("n=", nrow(orthologs.filtered), " orthologues"), 
       x = opt$name_ref, 
       y = opt$name_qry)+
  geom_hline(yintercept = c(0,fai_query$y.chrend/10^6), 
             linetype = "dotted", color = "gray50")+
  geom_vline(xintercept = c(0,fai_ref$x.chrend/10^6), 
             linetype = "dotted", color = "gray50")+
  theme(axis.line = element_line(linetype = "solid", size = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.title.x.bottom = element_text(size = 14, face = "italic"), 
        axis.title.y.left = element_text(size = 14, face = "italic"),
        axis.title.x.top = element_text(size = 14),
        axis.title.y.right = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"), 
        axis.text.x.bottom = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "right")

p

filename=paste("Oxford_plot", opt$output_filename, opt$locustag_ref, opt$locustag_query, sep = ".")
ggsave(paste(filename,"png",sep = "."), width = opt$plot_size, height = opt$plot_size, units = "in", dpi = 300)
ggsave(paste(filename,"tiff",sep = "."), width = opt$plot_size, height = opt$plot_size, units = "in", dpi = 300)

