base=/home/moritz/Data/people/0006_Erken_Genome_Project/1000_SAGs

usearch9 -sintax $base/1000_Data/2017-05-30_Erken\ univ\ 16S\ PCR\ screening.fasta -db ~/Data/utax/modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.sub.udb -tabbedout $base/1000_Data/2017-05-30_Erken_SAGs.sintax -strand both -sintax_cutoff 0.6

grep Bacteroidetes $base/1000_Data/2017-05-30_Erken_SAGs.sintax | grep -v Flavo | cut -f1 > $base/1000_Data/non_flavo_bacts.txt
grep Bacteroidetes $base/1000_Data/2017-05-30_Erken_SAGs.sintax | grep Flavo | grep 'd:Bacteria,p:Bacteroidetes$\|d:Bacteria$' | cut -f1 >> $base/1000_Data/non_flavo_bacts.txt
