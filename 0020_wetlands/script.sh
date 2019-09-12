base_path=$HOME/people/0020_wetlands
raws_folder=$base_path/000_data/020_reads/

python_hjelper=$HOME/Dropbox/Documents/Uppsala/Limno/0002_stump/helper.py
usearch=usearch10
usearch_path=$base_path/200_usearch/
nproc=11

mkdir $base_path
cd  $base_path

mkdir $base_path/000_data

rsync -aP moritz@rackham.uppmax.uu.se:/home/moritz/proj_folder/uppstore2017149/b2015084/INBOX/150824_M00485_0221_000000000-AG13K/ $base_path/000_data/020_reads/

mkdir $base_path/100_formated_data


samples=($(ls  $raws_folder | grep Sample | grep  WU | cut -f2 -d"_" | sort | uniq))

mkdir $base_path/100_formated_data/merged/

for s in ${samples[@]}
do
  echo $s
  fwd=`ls $raws_folder/Sample_$s | grep "_R1_" | grep -v fastqc`
  rev=`ls $raws_folder/Sample_$s | grep "_R2_" | grep -v fastqc`
  $usearch -fastq_mergepairs $raws_folder/Sample_${s}/$fwd -reverse $raws_folder/Sample_${s}/$rev -fastqout $base_path/100_formated_data/merged/${s}.fastq --relabel "$s:" -fastq_maxdiffs 10 -fastq_trunctail 5 -fastq_pctid 80  2>&1 | grep "% merged" | rev | cut -f2 -d" " | rev | sed "s/^/$s\t/" >> output_merging.txt
done
#
cat $base_path/100_formated_data/merged/*.fastq > $base_path/100_formated_data/all_reads.fastq
mkdir $usearch_path

$usearch -fastq_filter $base_path/100_formated_data/all_reads.fastq -fastq_minlen 250 -fastq_maxee  2 -fastaout $base_path/100_formated_data/all_qced_reads.fasta
$usearch -fastq_filter $base_path/100_formated_data/all_reads.fastq -fastaout $base_path/100_formated_data/all_reads.fasta

#$usearch --fastx_uniques $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --fastaout $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH

$usearch --fastx_uniques $base_path/test.mongo  -relabel Uniq_ -sizeout --minuniquesize 2 --fastaout $usearch_path/uniques_nosingletons.fa

$usearch -cluster_otus  $usearch_path/uniques_nosingletons.fa -otus $usearch_path/otus.fasta
$usearch -usearch_global $base_path/100_formated_data/all_qced_reads.fasta -db $usearch_path/otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/otu_table.txt
$usearch -sintax $usearch_path/otus.fasta -db ~/Data/tax_dbs/rdp_16s_v16_sp.udb -tabbedout $usearch_path/usearch.sintax -strand both


for s in ${samples[@]}
do
  echo $s
  $usearch -usearch_global $base_path/100_formated_data/merged/${s}.fastq -db $usearch_path/otus.fasta -strand both -id 0.97 2>&1 | grep "% matched" | rev | cut -f2 -d" " | rev | sed "s/^/$s\t/" >> mapping_rates.txt
done

### generate tables:

mkdir $base_path/300_tables

python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch_taxonomy.tax $base_path/300_tables/usearch_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise_taxonomy.tax $base_path/300_tables/unoise_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch.sintax $base_path/300_tables/usearch_sintax_silva 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise.sintax $base_path/300_tables/unoise_sintax_silva 0.8


mkdir $base_path/400_post_processing
