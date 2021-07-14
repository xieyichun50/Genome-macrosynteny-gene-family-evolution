#source /home/wyn/cgat-install/conda-install/bin/activate cgat-s
#cgat gtf2gtf -I ref_Adig_1.1_scaffolds.gtf.CDS  --method=filter --filter-method=longest-gene > ref_Adig_1.1_scaffolds.gtf.CDS.longest-gene.gtf
cat namelist | while read i;
do 
echo "$i"
/data/bin/gff2gtf.pl $i*.gff | grep -w CDS | sed 's/_id ".*"; transcript//' > $i.CDS;
cgat gtf2gtf -I $i.CDS  --method=filter --filter-method=longest-gene > $i.longest-gene.gtf
grep -w "gene_id" $i.longest-gene.gtf |sed 's/.*gene_id "//g;s/".*//' |sort -u > $i.longest-gene.gtf.id
seqtk subseq $i*.faa  $i.longest-gene.gtf.id > $i.longest-gene.proteins.fa
wc -l $i.longest-gene.gtf.id
wc -l $i.longest-gene.proteins.fa
done

