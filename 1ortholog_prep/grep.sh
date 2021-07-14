cat list | while read i;
do
grep "$i" Biomphalaria_straminea_hic.proteins.fa.filtered.fa.fai >> deg.list
done

cut -f1,2 deg.list | sed 's/\t/:1-/' > deg.list1

cat deg.list1 | while read i;
do
samtools faidx Biomphalaria_straminea_hic.proteins.fa.filtered.fa "$i" >> deg.nucl.fa
done
