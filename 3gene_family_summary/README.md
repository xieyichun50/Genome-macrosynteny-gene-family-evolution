# Summarise gene families
Procedure:
1. Run `1prepcafe.sh` and `2runcafe.sh`
2. Run `perl orthofinder_parser_genenum.tab.pl` in the same directory of `OrthoFinder/$date/Orthogroups/Orthogroups.GeneCount.tsv`
3. Gene family summary on tree node and tips + ortholog stats: `merge_tree_myriapod.R`
4. heatmap on shared orthologues `orthofinder_overlap_myriapod.R`

Files:
1. `SpeciesTree_rooted_node_labels.txt` in `Orthofinder/$date/Species_Tree` folder.
2. `SpeciesTreeAlignment.fa` in `OrthoFinder/$date/MultipleSequenceAlignments` folder
3. `OrthoFinder/$date/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv`
4. `OrthoFinder/$date/Orthogroups/Orthogroups.GeneCount.tsv`

Software:
1. [r8s](https://sourceforge.net/projects/r8s/files/r8s1.81.tar.gz)
2. [CAFE5](https://github.com/hahnlab/CAFE5)
3. CAFE python scripts: you may find from `cafe_pythonscript` folder
4. samtools
5. R/Rscript
6. R packages: 'ape', 'ggplot2', 'tidytree', 'ggtree', 'flextable', 'tidyr', 'dplyr', 'stringr', 'svglite', 'ggplotify', 'eoffice' 
