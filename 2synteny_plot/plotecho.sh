wd=/jelly_data/yichun/myriapod
ls | grep "ortholog" | grep "txt" | while read i;
do
xspecies=$(echo $i | sed 's/ortholog.//;s/__v__/\t/;s/.tsv.txt//' | cut -f1)
yspecies=$(echo $i | sed 's/ortholog.//;s/__v__/\t/;s/.tsv.txt//' | cut -f2)
echo "echo \"$xspecies and $yspecies\""
echo "Rscript /jelly_data/yichun/scripts/2synteny/orthofinder_synteny_plot_ortho_frag_significance_longgtf.R --input $i --xspecies-fai ${wd}/fai/${xspecies}.fa.fai --yspecies-fai ${wd}/fai/${yspecies}.fa.fai --xspecies-gtf ${wd}/gtf/formatted.${xspecies}.longest-gene.gtf --yspecies-gtf ${wd}/gtf/formatted.${yspecies}.longest-gene.gtf --xspecies-name ${xspecies} --yspecies-name ${yspecies} --min-orthologs-x 20 --min-orthologs-y 5 --plot-size 15"
done
