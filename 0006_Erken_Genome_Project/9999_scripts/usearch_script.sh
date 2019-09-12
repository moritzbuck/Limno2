base_path=$HOME/people/0006_Erken_Genome_Project/20000_mesocosm_amplicon_1/
raws_folder=/home/moritz/people/0006_Erken_Genome_Project/0000_rawdata/0200_mesocosm_amplicon_I/

python_hjelper=/home/moritz/temp/helper.py
usearch=usearch10
usearch_path=$base_path/2000_usearch/
nproc=19
db=~/dbs/rdp_16s_v16_sp.udb

source ~/people/0010_Pauline/arctic/bin/activate

cd  $base_path

#python ~/temp/fastq_deplex.py $raws_folder/Undetermined/Undetermined_S0_L001_R1_001.fastq.gz $raws_folder/Undetermined/Undetermined_S0_L001_R2_001.fastq.gz  /home/moritz/people/0024_haiyan/0000_data/haiyan_bcs.txt $raws_folder/deplexed/ $raws_folder/deplex_stats.txt False False

name=mesocosm_all

mkdir -p $base_path/1000_formated_data/1100_$name/
touch $base_path/1000_formated_data/output_merging.txt

samples=($(find  $raws_folder -name "*.fastq"  | grep _R1  |  rev | cut -f1 -d"/" | rev | cut -f2 -d_))
echo "sample" $'\t' "merged_reads" > $base_path/1000_formated_data/1100_$name/output_merging.txt

for s in ${samples[@]}
do
echo $s
  fwd=`ls $raws_folder/P11402_${s}_*_R1_001.fastq`
  rev=`ls $raws_folder/P11402_${s}_*_R2_001.fastq`
  $usearch -threads 16 -fastq_mergepairs $fwd -reverse $rev -fastqout $base_path/1000_formated_data/${s}.fastq --relabel "$s:" -fastq_maxdiffs 50 -fastq_trunctail 5  2>&1 | grep -E '% merged' | tail -n1  | sed "s/.* \([0-9.].*%\) merged/$s\t\1/"  >> $base_path/1000_formated_data/output_merging.txt
done
cat $base_path/1000_formated_data/*.fastq  > $base_path/1000_formated_data/all.fastq

mkdir -p $usearch_path/

vsearch -fastq_filter $base_path/1000_formated_data/all.fastq -fastq_minlen 400 -fastq_maxee  1 -fastaout $base_path/1000_formated_data/all.qced.fasta
vsearch -fastq_filter $base_path/1000_formated_data/all.fastq -fastaout $base_path/1000_formated_data/all.fasta
#sed -i 's/-/_/g' $base_path/1000_formated_data/all.fasta

vsearch --derep_fulllength  $base_path/1000_formated_data/all.qced.fasta -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/uniques_nosingletons.fa
$usearch -cluster_otus  $usearch_path/uniques_nosingletons.fa -otus $usearch_path/otus.fasta

vsearch -usearch_global $base_path/1000_formated_data/all.fasta -db $usearch_path/otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/otu_table.txt

$usearch -sintax $usearch_path/otus.fasta -db $db -tabbedout $usearch_path/taxonomy.sintax -strand both


#$usearch --fastx_uniques $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --fastaout $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH


### generate tables:

mkdir -p $base_path/3000_tables

python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/taxonomy.sintax $base_path/3000_tables/ 0.8


mkdir -p $base_path/4000_post_processing
