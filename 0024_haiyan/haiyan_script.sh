base_path=$HOME/people/0024_haiyan
raws_folder=$base_path/0000_data/0001_reads/180813_M00485_0446_000000000-BJDBP/

python_hjelper=$base_path/helper.py
usearch=usearch10
usearch_path=$base_path/2000_usearch/
nproc=11
db=~/dbs/rdp_16s_v16_sp.udb

source ~/people/0010_Pauline/arctic/bin/activate

cd  $base_path

python ~/temp/fastq_deplex.py $raws_folder/Undetermined/Undetermined_S0_L001_R1_001.fastq.gz $raws_folder/Undetermined/Undetermined_S0_L001_R2_001.fastq.gz  /home/moritz/people/0024_haiyan/0000_data/haiyan_bcs.txt $raws_folder/deplexed/ $raws_folder/deplex_stats.txt False False

name=haiyan

mkdir -p $base_path/1000_formated_data/1100_$name/
touch $base_path/1000_formated_data/1100_$name/output_merging.txt

samples=($(find  $raws_folder/deplexed/ -name "*.fastq"  | grep _fwd  | rev | cut -f2 -d_ | cut -f1 -d"/" | rev))
echo "sample" $'\t' "merged_reads" > $base_path/1000_formated_data/1100_$name/output_merging.txt

for s in ${samples[@]}
do
echo $s
  fwd=$raws_folder/deplexed/$s/${s}_fwd.fastq
  rev=$raws_folder/deplexed/$s/${s}_rev.fastq
  $usearch -fastq_mergepairs $fwd -reverse $rev -fastqout $base_path/1000_formated_data/1100_$name/${s}.fastq --relabel "$s:" -fastq_maxdiffs 50 -fastq_trunctail 5  2>&1 | grep -E '% merged' | tail -n1  | sed "s/.* \([0-9.].*%\) merged/$s\t\1/"  >> $base_path/1000_formated_data/1100_$name/output_merging.txt
done
cat $base_path/1000_formated_data/1100_$name/*.fastq  > $base_path/1000_formated_data/$name.fastq

mkdir -p $usearch_path/2100_$name

vsearch -fastq_filter $base_path/1000_formated_data/$name.fastq -fastq_minlen 400 -fastq_maxee  2 -fastaout $base_path/1000_formated_data/$name.qced.fasta
vsearch -fastq_filter $base_path/1000_formated_data/$name.fastq -fastaout $base_path/1000_formated_data/$name.fasta
sed -i 's/-/_/g' $base_path/1000_formated_data/$name.fasta

vsearch --derep_fulllength  $base_path/1000_formated_data/$name.qced.fasta -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/2100_$name/uniques_nosingletons.fa
$usearch -cluster_otus  $usearch_path/2100_$name/uniques_nosingletons.fa -otus $usearch_path/2100_$name/otus.fasta

vsearch --derep_fulllength  $base_path/1000_formated_data/$name.qced.fasta -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/2100_$name/uniques_w_singletons.fa
$usearch -cluster_otus  $usearch_path/2100_$name/uniques_w_singletons.fa -otus $usearch_path/2100_$name/otus_w_singletons.fasta


$usearch -usearch_global $base_path/1000_formated_data/$name.fasta -db $usearch_path/2100_$name/otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/2100_$name/otu_table.txt

$usearch -usearch_global $base_path/1000_formated_data/$name.fasta -db $usearch_path/2100_$name/otus_w_singletons.fasta -strand plus -id 0.97 -otutabout $usearch_path/2100_$name/otu_table_w_singletons.txt

$usearch -sintax $usearch_path/2100_$name/otus.fasta -db $db -tabbedout $usearch_path/2100_$name/taxonomy.sintax -strand both


#$usearch --fastx_uniques $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --fastaout $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH


### generate tables:

mkdir -p $base_path/3000_tables

python $python_hjelper make_taxo_table $usearch_path/2100_$name/otu_table.txt $usearch_path/2100_$name/taxonomy.sintax $base_path/3000_tables/3100_$name/ 0.8


mkdir -p $base_path/4000_post_processing
