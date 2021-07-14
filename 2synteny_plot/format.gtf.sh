cat /jelly_data/yichun/snail/longest_protO/specieslist | while read i;
do 
Rscript /jelly_data/yichun/scripts/2synteny/orthofinder_gtf_formatting.R -i $i.longest-gene.gtf 
done
