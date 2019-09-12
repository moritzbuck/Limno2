#!/usr/bin/bash

data_path=$HOME/DataBases/raws/100cultures
data_path=$HOME/Data/100cultures
usearch=usearch10.0.240_i86linux32
nproc=11
bbmerge=/home/moritz/src/bbmap/bbmerge.sh


gunzip `find $data_path -name "*.gz"`

for d in `find $data_path -type d -name "Sample_*"`
do
    sample=`basename $d`
    sample=${sample##Sample_}
    sample=`echo $sample | sed 's/Erken/erken/'`
    dd=`echo $sample | sed 's/-/_/'`
    echo $sample
#    pandaseq -F -T $nproc -f $d/*_R1_*.fastq -r $d/*_R2_*.fastq -w $d/${sample}_merged_panda.fastq  2> /dev/null
    $usearch -fastq_mergepairs $d/*_R1_*.fastq -reverse $d/*_R2_*.fastq -fastqout $d/${sample}.fastq --relabel "$sample:" -fastq_maxdiffs 50  2>&1 | grep "% merged" | rev | cut -f2 -d" " | rev | sed "s/^/$sample\t/" >> output_merging.txt
    $usearch -fastq_filter $d/${sample}.fastq -fastaout $d/${sample}_unfil.fasta > /dev/null 2> /dev/null
   $usearch -fastq_filter $d/${sample}.fastq -fastq_maxee 1 -fastq_minlen 150 -fastaout $d/${sample}.fasta 2>&1 | grep "% passed" | rev | cut -f2 -d" " | rev | sed "s/^/$sample\t/" >> output_qc.txt
done

cat `find $data_path -name "00*.fasta" | grep -v _unfil | grep -v _PLUS-C` >  all_negatives.fasta

cat `find $data_path -name "[Ee]rken*.fasta" | grep -v _unfil` >  all_erken.fasta
cat `find $data_path -name "[Ee]rken*.fasta" | grep _unfil` >  all_erken_unfil.fasta

sed -i 's/-/_/' all_erken.fasta
sed -i 's/-/_/' all_negatives.fasta

cat all_negatives.fasta >> all_erken.fasta

$usearch -fastx_uniques  all_erken.fasta -relabel Uniq -sizeout --minuniquesize 2 --fastaout  erken_uniques_nosingletons.fa

#  52 OTUs 49.8 %  OTUs in NTC


$usearch -cluster_otus erken_uniques_nosingletons.fa -otus erken_otus_97_nosingletons.fa
$usearch -usearch_global all_erken.fasta -db erken_otus_97_nosingletons.fa -strand plus -id 0.97 -otutabout erken_otu_table_97_nosingletons.txt

$usearch -sintax erken_otus_97_nosingletons.fa -db ~/Data/tax_dbs/rdp_16s_v16_sp.v9.udb -tabbedout erken_97_singletons.rdp.sintax -strand both -sintax_cutoff 0.8
# sed -i 's/\t/;/g' erken_otus_97_nosingletons.rdp.tax sed -i 's/\t/;/g' erken_otus_97_nosingletons.rdp.tax

# 144 OTUs 51.1% 42 OTUs in NTC, 74 in 00


$usearch -unoise2 erken_uniques_nosingletons.fa -fastaout erken_otus100_nosingletons.fa -minampsize 4 -ampout erken_otus100_nosingletons.allamps.fa
$usearch -usearch_global all_erken.fasta -db erken_otus100_nosingletons.fa -strand plus -id 0.97 -otutabout erken_otu_table_97_unoise_nosingletons.txt
$usearch -utax erken_otus100_nosingletons.fa -db ~/DataBases/utax/refdb_16.udb -utaxout erken_100.tax -strand both
$usearch -sintax erken_otus100_nosingletons.fa -db ~/DataBases/utax/modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta -tabbedout erken_100.sintax -strand both -sintax_cutoff 0.8

# 40 OTUs 49.1 %


# JAMTLAND

cat `find $data_path -name "Jamtla*.fasta" | grep -v _unfil` >  all_jamt.fasta
cat `find $data_path -name "Ana*.fasta" | grep -v _unfil` >> all_jamt.fasta
sed -i 's/-/_/' all_jamt.fasta
cat all_negatives.fasta >> all_jamt.fasta

$usearch --derep_fulllength all_jamt.fasta -relabel Uniq -sizeout --minuniquesize 3 --output  jamt_uniques_nodoublets.fa
$usearch --derep_fulllength all_jamt.fasta -relabel Uniq -sizeout --minuniquesize 2 --output  jamt_uniques_nosingletons.fa

$usearch -fastx_uniques  all_jamt.fasta -relabel Uniq -sizeout --minuniquesize 2 --fastaout  jamt_uniques_nosingletons.fa
$usearch -cluster_otus jamt_uniques_nosingletons.fa -otus jamt_otus_97_nosingletons.fa
$usearch -usearch_global all_jamt.fasta -db jamt_otus_97_nosingletons.fa -strand plus -id 0.97 -otutabout jamt_otu_table_97_nosingletons.txt

usearch9.2.64_i86linux32 -sintax jamt_otus_97_nosingletons.fa -db ~/Data/tax_dbs/rdp_16s_v16_sp.v9.udb -tabbedout jamt_97_singletons.rdp.sintax -strand both -sintax_cutoff 0.8
