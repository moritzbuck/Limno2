base_path=/home/moritz/people/0007_stefans_cores
raws_folder=$base_path/000_data/010_reads/
python_hjelper=$HOME/Dropbox/Documents/Uppsala/Limno/0007_stefans_cores/helper.py
usearch=usearch9
usearch_path=$base_path/200_usearch/
nproc=11

cd  $base_path
samples=($(ls $base_path/100_formated_data/deplexed | sed 's/_fwd.fastq//' | sed 's/_rev.fastq//' | sort | uniq))

mkdir $base_path/100_formated_data


for f in `ls $base_path/000_data/020_jsons/*.json | grep DE21`
do
  python $python_hjelper deplex_reads $f $raws_folder $base_path/100_formated_data/deplexed/ 7
done

find $base_path/100_formated_data/deplexed -type f -name "* *.fastq" -exec rename "s/\s/_/g" {} \;
samples=($(ls $base_path/100_formated_data/deplexed | sed 's/_fwd.fastq//' | sed 's/_rev.fastq//' | sort | uniq))

mkdir $base_path/100_formated_data/merged/
pp=$base_path/100_formated_data/deplexed

for s in ${samples[@]}
do
  echo $s
  $usearch -fastq_mergepairs $pp/${s}_fwd.fastq -reverse $pp/${s}_rev.fastq -fastqout $base_path/100_formated_data/merged/${s}.fastq --relabel "$s:" -fastq_maxdiffs 15
done > temp.out 2> temp.err

cat $base_path/100_formated_data/merged/*.fastq > $base_path/100_formated_data/all_reads.fastq
mkdir $usearch_path


vsearch -fastq_filter $base_path/100_formated_data/all_reads.fastq -fastq_maxee 1 -fastaout $base_path/100_formated_data/all_qced_reads.fasta
vsearch -fastq_filter $base_path/100_formated_data/all_reads.fastq -fastaout $base_path/100_formated_data/all_reads.fasta
vsearch --derep_fulllength $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH

$usearch -cluster_otus  $usearch_path/uniques_nosingletons.fa -otus $usearch_path/otus.fasta
vsearch -usearch_global $base_path/100_formated_data/all_reads.fasta -db $usearch_path/otus2.fasta -strand plus -id 0.97 -otutabout $usearch_path/otu_table2.txt >
$usearch -utax $usearch_path/otus.fasta -db ~/DataBases/utax/refdb_16.udb -utaxout $usearch_path/usearch_taxonomy.tax -strand both
$usearch -sintax $usearch_path/otus.fasta -db ~/DataBases/utax/modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta -tabbedout $usearch_path/usearch.sintax -strand both -sintax_cutoff 0.8
for s in ${samples[@]}; do
echo $s `usearch9 -usearch_global /home/moritz/people/0007_stefans_cores/100_formated_data/merged/${s}.fastq -db $usearch_path/otus.fasta -strand plus -id 0.97 2>&1 > /dev/null | grep matched | rev | cut -d',' -f1 | rev`;
done >  $usearch_path/usearch_mapping_rates_persample.txt

### UNOISE

$usearch -unoise2 $usearch_path/uniques_nosingletons.fa -fastaout $usearch_path/unoise_otus.fasta -minampsize 4
vsearch -usearch_global $base_path/100_formated_data/all_reads.fasta -db $usearch_path/unoise_otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/unoise_otu_table.txt
$usearch -utax $usearch_path/unoise_otus.fasta -db ~/DataBases/utax/refdb_16.udb -utaxout $usearch_path/unoise_taxonomy.tax -strand both
$usearch -sintax $usearch_path/unoise_otus.fasta -db ~/DataBases/utax/modified_SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta -tabbedout $usearch_path/unoise.sintax -strand both -sintax_cutoff 0.8

### generate tables:

mkdir $base_path/300_tables

python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch_taxonomy.tax $base_path/300_tables/usearch_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise_taxonomy.tax $base_path/300_tables/unoise_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch.sintax $base_path/300_tables/usearch_sintax_silva 0.8


mkdir $base_path/400_post_processing


sh.raxmlHPC_AVX('-m', "PROTGAMMALGF", "-T", threads-2 , '-p', sede, '-s', alignment, '-n', 'tree', '-w', rax_path, '-f', 'a', '-x', 1, '-N', 'autoMR')
