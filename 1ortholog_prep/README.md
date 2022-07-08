# Prepare longest protein sequence fasta file and ortholog calling.

Procedure:
1. run `sh filter_longest_genes_from_gff*.sh` to get the longest protein. Prepare a `namelist3` file with species name, or do `cut -f1 links.txt > namelist3`.
2. run `sh orthofinder.sh` to perform ortholog analysis

Files:
1. genome.gff/gff3: gene annotation file in gff or gff3 format
2. protein.fasta: all protein sequence in fasta format

Software:
1. samtools
2. [seqtk](https://github.com/lh3/seqtk)
3. [cgat](https://github.com/cgat-developers/cgat-apps)
4. [orthofinder](https://github.com/davidemms/OrthoFinder)

For some gff files, the auto script might not work. Prepare the `*.CDS` file in the following format and continue with `cgat` line.
```
NC_007419.2	Gnomon	CDS	3419327	3419833	0	+	0	gene_id "XP_008191933.1"; 
NC_007419.2	Gnomon	CDS	3421002	3421537	0	+	0	gene_id "XP_008191933.1"; 
NC_007419.2	Gnomon	CDS	3421594	3421693	0	+	1	gene_id "XP_008191933.1"; 
NC_007419.2	Gnomon	CDS	3421746	3421812	0	+	0	gene_id "XP_008191933.1"; 
NC_007419.2	Gnomon	CDS	3421863	3421942	0	+	2	gene_id "XP_008191933.1"; 
NC_007422.5	Gnomon	CDS	8428903	8428975	0	+	0	gene_id "XP_008195602.1"; 
NC_007422.5	Gnomon	CDS	8442559	8442600	0	+	2	gene_id "XP_008195602.1"; 
NC_007422.5	Gnomon	CDS	8475212	8475378	0	+	2	gene_id "XP_008195602.1"; 
NC_007422.5	Gnomon	CDS	8477468	8477569	0	+	0	gene_id "XP_008195602.1"; 
NC_007422.5	Gnomon	CDS	8480494	8480602	0	+	0	gene_id "XP_008195602.1"; 
```
