cat genefile | while read i;
do 
sed 's/ /_/g;s/(/_/g;s/)/_/g' $i.fas > $i.align.fa
echo "iqtree -s $i.align.fa --seqtype AA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000"
iqtree -s $i.align.fa --seqtype AA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000 
done
