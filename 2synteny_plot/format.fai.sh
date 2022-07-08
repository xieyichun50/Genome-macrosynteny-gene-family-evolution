
cat specieslist | while read i;
do 
samtools faidx $i.fa
done
