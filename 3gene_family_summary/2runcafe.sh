#run either of the following regarding species

#11 species, 9 myriapod + 2 horseshoe crab
#python2 /tools/CAFE5/python_scripts/cafetutorial_prep_r8s.py -i SpeciesTree_rooted_node_labels.txt -o r8s_ctl_file.txt -s 671294 -p 'Carcinoscorpius_rotundicauda,Rhysida_immarginata' -c '601'

#24 species, no human
#python2 /tools/CAFE5/python_scripts/cafetutorial_prep_r8s.py -i SpeciesTree_rooted_node_labels.txt -o r8s_ctl_file.txt -s 156815 -p 'Caenorhabditis_elegans,Branchiostoma_floridae' -c '736'

#23 species, no human no Branchiostoma_floridae
#python2 /tools/CAFE5/python_scripts/cafetutorial_prep_r8s.py -i SpeciesTree_rooted_node_labels.txt -o r8s_ctl_file.txt -s 166951 -p 'Caenorhabditis_elegans,Drosophila_melanogaster' -c '728'

#25 species, with human and Branchiostoma_floridae
#python2 /tools/CAFE5/python_scripts/cafetutorial_prep_r8s.py -i SpeciesTree_rooted_node_labels.txt -o r8s_ctl_file.txt -s 138195 -p 'Caenorhabditis_elegans,Homo_sapiens' -c '736'

/tools/r8s/r8s1.81/src/r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > r8s.ultrametric.tre
cat r8s.ultrametric.tre

##run cafe
echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda1 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda1 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda2 -k 2 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda2 -k 2 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda3 -k 3 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda3 -k 3 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda4 -k 4 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda4 -k 4 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda5 -k 5 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda5 -k 5 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda6 -k 6 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda6 -k 6 -c 180

#echo "/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda1 -c 180"
#/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda1 -c 180

#echo "/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda2 -k 2 -c 180"
#/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda2 -k 2 -c 180

#echo "/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda3 -k 3 -c 180"
#/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda3 -k 3 -c 180

#echo "/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda4 -k 4 -c 180"
#/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda4 -k 4 -c 180

#echo "/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda5 -k 5 -c 180"
#/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda5 -k 5 -c 180

#echo "/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda6 -k 6 -c 180"
#/tools/CAFE5/bin/cafe5 -i cafe.input.tsv -t r8s.ultrametric.tre -o r8s_lambda6 -k 6 -c 180
