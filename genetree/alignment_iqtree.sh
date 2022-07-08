cat genefile | while read i;
do
echo "Aligning $i"
mafft --thread 38 --auto $i.fa > $i.align.fa
echo "iqtree -s $i.align.fa --seqtype AA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000"
iqtree -s $i.align.fa --seqtype AA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000 
done
