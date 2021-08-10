find . | xargs ls | grep "Orthologues" | grep "tsv" | while read i;
do
mv $i .
done

ls | grep "tsv" | while read i;
do
echo "Rscript /jelly_data/yichun/scripts/2synteny/orthofinder_orthologues_formating.R -i $i "
Rscript /jelly_data/yichun/scripts/2synteny/orthofinder_orthologues_formating.R -i $i
done
