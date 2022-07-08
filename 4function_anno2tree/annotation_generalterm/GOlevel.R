suppressPackageStartupMessages(library(annotate))
suppressPackageStartupMessages(library(GO.db))

getAllBPChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

getAllMFChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOMFCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

getAllCCChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOCCCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

##MF
level3_terms.MF <- getAllMFChildren("GO:0003674")
level4_terms.MF <- getAllMFChildren(level3_terms.MF)
level5_terms.MF <- getAllMFChildren(level4_terms.MF)
level6_terms.MF <- getAllMFChildren(level5_terms.MF)
level7_terms.MF <- getAllMFChildren(level6_terms.MF)

level3_terms.MF<-as.data.frame(level3_terms.MF)
names(level3_terms.MF)[1]="goClass"

level4_terms.MF<-as.data.frame(level4_terms.MF)
names(level4_terms.MF)[1]="goClass"

level5_terms.MF<-as.data.frame(level5_terms.MF)
names(level5_terms.MF)[1]="goClass"

level6_terms.MF<-as.data.frame(level6_terms.MF)
names(level6_terms.MF)[1]="goClass"

level7_terms.MF<-as.data.frame(level7_terms.MF)
names(level7_terms.MF)[1]="goClass"

##BP
level3_terms.BP <- getAllBPChildren("GO:0002376")
level4_terms.BP <- getAllBPChildren(level3_terms.BP)
level5_terms.BP <- getAllBPChildren(level4_terms.BP)
level6_terms.BP <- getAllBPChildren(level5_terms.BP)
level7_terms.BP <- getAllBPChildren(level6_terms.BP)

level3_terms.BP<-as.data.frame(level3_terms.BP)
names(level3_terms.BP)[1]="goClass"

level4_terms.BP<-as.data.frame(level4_terms.BP)
names(level4_terms.BP)[1]="goClass"

level5_terms.BP<-as.data.frame(level5_terms.BP)
names(level5_terms.BP)[1]="goClass"

level6_terms.BP<-as.data.frame(level6_terms.BP)
names(level6_terms.BP)[1]="goClass"

level7_terms.BP<-as.data.frame(level7_terms.BP)
names(level7_terms.BP)[1]="goClass"

##CC
level3_terms.CC <- getAllCCChildren("GO:0005575")
level4_terms.CC <- getAllCCChildren(level3_terms.CC)
level5_terms.CC <- getAllCCChildren(level4_terms.CC)
level6_terms.CC <- getAllCCChildren(level5_terms.CC)
level7_terms.CC <- getAllCCChildren(level6_terms.CC)

level3_terms.CC<-as.data.frame(level3_terms.CC)
names(level3_terms.CC)[1]="goClass"

level4_terms.CC<-as.data.frame(level4_terms.CC)
names(level4_terms.CC)[1]="goClass"

level5_terms.CC<-as.data.frame(level5_terms.CC)
names(level5_terms.CC)[1]="goClass"

level6_terms.CC<-as.data.frame(level6_terms.CC)
names(level6_terms.CC)[1]="goClass"

level7_terms.CC<-as.data.frame(level7_terms.CC)
names(level7_terms.CC)[1]="goClass"

##rbind
level3_terms<-rbind(level3_terms.MF, level3_terms.BP, level3_terms.CC)
level3_terms$level="L3"

level4_terms<-rbind(level4_terms.MF, level4_terms.BP, level4_terms.CC)
level4_terms$level="L4"

level5_terms<-rbind(level5_terms.MF, level5_terms.BP, level5_terms.CC)
level5_terms$level="L5"

level6_terms<-rbind(level6_terms.MF, level6_terms.BP, level6_terms.CC)
level6_terms$level="L6"

level7_terms<-rbind(level7_terms.MF, level7_terms.BP, level7_terms.CC)
level7_terms$level="L7"

write.table(level3_terms, file = "D:/3enrichment/GOlevel3.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

write.table(level4_terms, file = "D:/3enrichment/GOlevel4.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

write.table(level5_terms, file = "D:/3enrichment/GOlevel5.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

write.table(level6_terms, file = "D:/3enrichment/GOlevel6.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

write.table(level7_terms, file = "D:/3enrichment/GOlevel7.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)

result %<>% subset(ID %in% terms)

