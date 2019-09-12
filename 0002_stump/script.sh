base_path=/home/moritzbuck/people/0002_stump
raws_folder=$base_path/000_data/020_reads/
hgca_folder=~/uppmax/people/jingjing
python_hjelper=$HOME/people/0002_stump/888_scripts/helper.py
usearch=usearch9
usearch_path=$base_path/200_usearch/
nproc=18

mkdir $base_path
cd  $base_path

mkdir $base_path/000_data

cp -r ~/uppmax/temp/jsons/ $base_path/000_data/010_jsons


mkdir $base_path/100_formated_data

mkdir 00A_hgcA_amplicon
cp $hgca_folder/hgca_amplicon/otu_table_fwd.csv 00A_hgcA_amplicon
cp $hgca_folder/hgca_amplicon/centers_fwd.fasta 00A_hgcA_amplicon
cp $hgca_folder/classification.txt 00A_hgcA_amplicon

python $python_hjelper prefilter_hgca_table_and_stuff

mkdir 00B_muscle

cat 00A_hgcA_amplicon/centers_stump.faa 000_data/HgcA_db.faa > 00B_muscle/all_hgcas.faa
cat 00A_hgcA_amplicon/centers_stump_lt1000.fasta 000_data/HgcA_db.faa > 00B_muscle/less_hgcas.faa
cat 00A_hgcA_amplicon/centers_hwu_lt1000.fasta 000_data/HgcA_db.faa > 00B_muscle/less_hwu.faa

muscle -in 00B_muscle/all_hgcas.faa -out 00B_muscle/HgcAs.aligned.faa
muscle -in 00B_muscle/less_hgcas.faa -out 00B_muscle/Less_HgcAs.aligned.faa
muscle -in 00B_muscle/less_hwu.faa -out 00B_muscle/Less_hwu.aligned.faa

mkdir 00C_raxml
mkdir 00D_less_raxml
mkdir 00E_hwu_raxml


raxmlHPC-PTHREADS -m PROTGAMMAIWAGX -T 10 -p 23 -s $base_path/00B_muscle/HgcAs.aligned.tunked.faa -n tree -w $base_path/00C_raxml -f a -x 1 -N autoMR
raxmlHPC-PTHREADS -m PROTGAMMAIWAGX -T 10 -p 23 -s $base_path/00B_muscle/Less_HgcAs.aligned.tunked.faa -n tree -w $base_path/00D_less_raxml -f a -x 1 -N autoMR
raxmlHPC-PTHREADS -m PROTGAMMAIWAGX -T 10 -p 23 -s $base_path/00B_muscle/Less_hwu.aligned.tunked.faa -n tree -w $base_path/00E_hwu_raxml -f a -x 1 -N autoMR


mkdir $base_path/100_formated_data/deplexed

for f in `ls $base_path/000_data/010_jsons/*.json`
do
  python $python_hjelper deplex_reads $f $raws_folder $base_path/100_formated_data/deplexed/ 7
done

samples=($(ls $base_path/100_formated_data/deplexed | cut -f1 -d"_" | sort | uniq))

mkdir $base_path/100_formated_data/merged/
pp=$base_path/100_formated_data/deplexed

for s in ${samples[@]}
do
  $usearch -fastq_mergepairs $pp/${s}_fwd.fastq -reverse $pp/${s}_rev.fastq -fastqout $base_path/100_formated_data/merged/${s}.fastq --relabel "$s:"
done

cat $base_path/100_formated_data/merged/stubb* > $base_path/100_formated_data/all_reads.fastq
mkdir $usearch_pathmkdir $usearch_path

$usearch -fastq_filter $base_path/100_formated_data/all_reads.fastq -fastq_maxee 1 -fastaout $base_path/100_formated_data/all_qced_reads.fasta
vsearch -fastq_filter $base_path/100_formated_data/all_reads.fastq -fastaout $base_path/100_formated_data/all_reads.fasta
vsearch --derep_fulllength $base_path/100_formated_data/all_qced_reads.fasta  -relabel Uniq_ -sizeout --minuniquesize 2 --output $usearch_path/uniques_nosingletons.fa

### STANDARD USEARCH

$usearch -cluster_otus  $usearch_path/uniques_nosingletons.fa -otus $usearch_path/otus.fasta
$usearch -usearch_global $base_path/100_formated_data/all_reads.fastq -db $usearch_path/otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/otu_table.txt
#$usearch -utax $usearch_path/otus.fasta -db ~/Data/utax/refdb_16.udb -utaxout $usearch_path/usearch_taxonomy.tax -strand both
$usearch -sintax $usearch_path/otus.fasta -db ~/data/dbs/SILVA_132_LSURef_SSURef_Nr99.fasta  -tabbedout $usearch_path/usearch.sintax -strand both -sintax_cutoff 0.8

### UNOISE

$usearch -unoise2 $usearch_path/uniques_nosingletons.fa -fastaout $usearch_path/unoise_otus.fasta -minampsize 4
vsearch -usearch_global $base_path/100_formated_data/all_reads.fasta -db $usearch_path/unoise_otus.fasta -strand plus -id 0.97 -otutabout $usearch_path/unoise_otu_table.txt
#$usearch -utax $usearch_path/unoise_otus.fasta -db ~/Data/utax/refdb_16.udb -utaxout $usearch_path/unoise_taxonomy.tax -strand both
$usearch -sintax $usearch_path/unoise_otus.fasta -db ~/data/dbs/SILVA_132_LSURef_SSURef_Nr99.fasta  -tabbedout $usearch_path/unoise.sintax -strand both -sintax_cutoff 0.8

### generate tables:

mkdir $base_path/300_tables

python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch_taxonomy.tax $base_path/300_tables/usearch_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise_taxonomy.tax $base_path/300_tables/unoise_utax_rdp_train 0.8
python $python_hjelper make_taxo_table $usearch_path/otu_table.txt $usearch_path/usearch.sintax $base_path/300_tables/usearch_sintax_silva 0.8
python $python_hjelper make_taxo_table $usearch_path/unoise_otu_table.txt $usearch_path/unoise.sintax $base_path/300_tables/unoise_sintax_silva 0.8


mkdir $base_path/400_post_processing
