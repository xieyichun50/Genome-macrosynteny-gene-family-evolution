ls | while read i ; 
do 
echo $i ; 
Rscript /jelly_data/yichun/scripts/2synteny/orthofinder_orthologues_formating.R -i $i ; 
done