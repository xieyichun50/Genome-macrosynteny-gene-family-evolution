##Read overlap count table
Orthofinder_dir="F:/Snail/longest_protO/Orthofinder/Results_May27/"
species.overlaps<-read.table(paste(Orthofinder_dir, 
                                   "Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv", 
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
nspecies=13
species.order<-subset(Species_tree_data, 
                      node <= nspecies,
                      select = c("y","label"))
species.order<-species.order[order(species.order$y),]$label
species.order
number.orthogroups=

heatmap.overlap.orthogroups<-ggplot(species.overlaps.1v1,
                                    aes(x=X, y=Y))+
  geom_tile(aes(fill = No.overlap.orthogroups/number.orthogroups))+
  geom_text(aes(label = No.overlap.orthogroups), color = "black", size = 4) +
  scale_fill_gradient(low = "white",high = "mediumpurple4"， limits =c(0,0.5))+
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
ggsave("heatmap.overlap.orthogroups.png", width = 12, height = 12, units = "in", dpi = 300)
ggsave("heatmap.overlap.orthogroups.tiff", width = 12, height = 12, units = "in", dpi = 300)
