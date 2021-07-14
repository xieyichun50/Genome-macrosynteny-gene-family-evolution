#source /home/wyn/cgat-install/conda-install/bin/activate cgat-s
#cgat gtf2gtf -I ref_Adig_1.1_scaffolds.gtf.CDS  --method=filter --filter-method=longest-gene > ref_Adig_1.1_scaffolds.gtf.CDS.longest-gene.gtf
cat namelist1 | while read i;
do 
#/tools/gffread/gffread $i*.gff3 -T -o $i.gtf 
echo "$i"
grep -w "CDS" $i*.gff | sed 's/_id ".*"://' | sed 's/transcript_id/\t/' | cut -f1,2,3,4,5,6,7,8,9  > $i.CDS;
cgat gtf2gtf -I $i.CDS  --method=filter --filter-method=longest-gene > $i.longest-gene.gtf
grep -w "gene_id" $i.longest-gene.gtf |sed 's/.*gene_id "//g;s/".*//' |sort -u > $i.longest-gene.gtf.id
seqtk subseq $i*.faa  $i.longest-gene.gtf.id > $i.longest-gene.proteins.fa
wc -l $i.longest-gene.gtf.id
wc -l $i.longest-gene.proteins.fa
done

