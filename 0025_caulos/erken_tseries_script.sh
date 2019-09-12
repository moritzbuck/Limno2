base_path=$HOME/people/0025_caulos
raws_folder=$base_path/0000_data/0100_erken_tseries/0110_rawdata

python_hjelper=$base_path/helper.py
usearch=usearch10
usearch_path=$base_path/2000_usearch/
nproc=11
db=~/dbs/rdp_16s_v16_sp.udb

mkdir -p $raws_folder

cd  $base_path

echo "/repository/user/main/public/root = '$raws_folder'" > $HOME/.ncbi/user-settings.mkfg
for f in `cat $raws_folder/SRR_Acc_List.txt `; do  prefetch -v $f ; done
fastq-dump --outdir $raws_folder/ --split-files $raws_folder/sra/*

name=erken_tseries
mkdir -p $base_path/1000_formated_data/1100_$name/
touch $base_path/1000_formated_data/1100_$name/output_merging.txt
samples=($(ls  $raws_folder | grep fastq | grep _1 | cut -f1 -d"_" | cut -f3- -d"-"))
echo "sample" $'\t' "merged_reads" > $base_path/1000_formated_data/1100_$name/output_merging.txt
for s in ${samples[@]}
do
echo $s
  fwd=$raws_folder/${s}_1.fastq
  rev=$raws_folder/${s}_2.fastq
  $usearch -fastq_mergepairs $fwd -reverse $rev -fastqout $base_path/1000_formated_data/1100_$name/${s}.fastq --relabel "$s:" -fastq_maxdiffs 50 -fastq_trunctail 5  2>&1 | grep -E '% merged' | tail -n1  | sed "s/.* \([0-9.].*%\) merged/$s\t\1/"  >> $base_path/1000_formated_data/1100_$name/output_merging.txt
done
cat $base_path/1000_formated_data/1100_$name/*.fastq  > $base_path/1000_formated_data/$name.fastq


mkdir -p $usearch_path/2100_$name

vsearch -fastq_filter $base_path/1000_formated_data/$name.fastq -fastq_minlen 400 -fastq_maxee  2 -fastaout $base_path/1000_formated_data/$name.qced.fasta
vsearch -fastq_filter $base_path/1000_formated_data/$name.fastq -fastaout $base_path/1000_formated_data/$name.fasta

vsearch --derep_fulllength  $base_path/1000_formated_data/$name.qced.fasta -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/2100_$name/uniques_nosingletons.fa
$usearch -cluster_otus  $usearch_path/2100_$name/uniques_nosingletons.fa -otus $usearch_path/2100_$name/otus.fasta

vsearch -usearch_global $base_path/1000_formated_data/$name.fasta -db $usearch_path/2100_$name/otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/2100_$name/otu_table.txt
$usearch -sintax $usearch_path/2100_$name/otus.fasta -db $db -tabbedout $usearch_path/2100_$name/taxonomy.sintax -strand both


#$usearch --fastx_uniques $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --fastaout $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH


### generate tables:

mkdir -p $base_path/3000_tables

python $python_hjelper make_taxo_table $usearch_path/2100_$name/otu_table.txt $usearch_path/2100_$name/taxonomy.sintax $base_path/3000_tables/3100_$name/ 0.8


mkdir $base_path/4000_post_processing
