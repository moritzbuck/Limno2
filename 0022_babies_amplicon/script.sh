base_path=$HOME/people/0022_babies_amplicon
raws_folder=$base_path/0000_data/0001_reads

python_hjelper=$base_path/helper.py
usearch=usearch10
usearch_path=$base_path/2000_usearch/
nproc=19
db=~/dbs/rdp_16s_v16_sp.udb


mkdir $base_path
cd  $base_path

ln -s $HOME/proj_folder/uppstore2017149/INBOX/RD-1801/180608_M00485_0430_000000000-BPM5J $raws_folder

unpigz $raws_folder/*.gz

mkdir $base_path/1000_formated_data
mkdir $base_path/1000_formated_data/Sanne
mkdir $base_path/1000_formated_data/Vaso
mkdir $base_path/1000_formated_data/Mayra

for name in `ls  $base_path/1000_formated_data | grep -v fast | grep Vaso`
do
  touch $base_path/1000_formated_data/$name/output_merging.txt
  samples=($(ls  $raws_folder | grep fastq | grep _R1_ | cut -f1 -d"_" | cut -f3- -d"-" | grep $name))
  echo "sample" $'\t' "merged_reads" > $base_path/1000_formated_data/$name/output_merging.txt
  for s in ${samples[@]}
  do
  echo $s
  fwd=`ls $raws_folder/RD*-$s\_* | grep "_R1_"`
  rev=`ls $raws_folder/RD*-$s\_* | grep "_R2_"`
  ss=`echo $s | tr "-" "_"`
  $usearch -fastq_mergepairs $fwd -reverse $rev -fastqout $base_path/1000_formated_data/$name/${s}.fastq --relabel "$ss:" -fastq_maxdiffs 5 -fastq_trunctail 5  2>&1 | grep -E '% merged' | tail -n1  | sed "s/.* \([0-9.].*%\) merged/$s\t\1/"  >> $base_path/1000_formated_data/$name/output_merging.txt
  done

  cat $base_path/1000_formated_data/$name/$name*.fastq  > $base_path/1000_formated_data/$name.fastq
done

name=Sanne
samples=($( ls $base_path/0000_data/0002_extra_Sanne  | grep fastq | grep _R1 | cut -f1-3 -d"_" ))
echo "sample" $'\t' "merged_reads" > $base_path/1000_formated_data/$name/output_merging.txt
for s in ${samples[@]}
do
echo $s
fwd=$base_path/0000_data/0002_extra_Sanne/${s}_R1.fastq
rev=$base_path/0000_data/0002_extra_Sanne/${s}_R2.fastq
$usearch -fastq_mergepairs $fwd -reverse $rev -fastqout $base_path/1000_formated_data/$name/${s}.fastq --relabel "$s:" -fastq_maxdiffs 5 -fastq_trunctail 5  2>&1 | grep -E '% merged' | tail -n1  | sed "s/.* \([0-9.].*%\) merged/$s\t\1/"  >> $base_path/1000_formated_data/$name/output_merging.txt
done
cat $base_path/1000_formated_data/$name/$name*.fastq  > $base_path/1000_formated_data/$name.fastq


mkdir $usearch_path

for f in `ls $base_path/1000_formated_data/*.fastq | grep Sanne`
do
name=`basename ${f%%.fastq}`
$usearch -fastq_filter $f -fastq_minlen 400 -fastq_maxee  1 -fastaout ${f%%.fastq}.qced.fasta
$usearch -fastq_filter $f -fastaout ${f%%.fastq}.fasta
mkdir $usearch_path/$name

vsearch --derep_fulllength  ${f%%.fastq}.qced.fasta -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/$name/uniques_nosingletons.fa
#$usearch -cluster_otus  $usearch_path/$name/uniques_nosingletons.fa -otus $usearch_path/$name/otus.fasta
vsearch --cluster_size  $usearch_path/$name/uniques_nosingletons.fa -id 0.97 --centroids $usearch_path/$name/otus.fasta

$usearch -usearch_global  $f -db $usearch_path/$name/otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/$name/otu_table.txt
$usearch -sintax $usearch_path/$name/otus.fasta -db $db -tabbedout $usearch_path/$name/taxonomy.sintax -strand both

done

#$usearch --fastx_uniques $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --fastaout $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH



for s in `ls 1000_formated_data/Sanne/*.fastq`
do
  name=Sanne
  $usearch -threads 20 -usearch_global $s -db $usearch_path/Sanne/otus.fasta -strand both -id 0.97 2>&1 | grep "% matched" | sed "s/.*Searching \(.*\).fastq, \([0-9.]*\)%.*/\1\t$name\t\2/"
done

### generate tables:

mkdir $base_path/3000_tables

python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch_taxonomy.tax $base_path/300_tables/usearch_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise_taxonomy.tax $base_path/300_tables/unoise_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch.sintax $base_path/300_tables/usearch_sintax_silva 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise.sintax $base_path/300_tables/unoise_sintax_silva 0.8


mkdir $base_path/400_post_processing
