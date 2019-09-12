



MY_PATH=/proj/g2017026/2017_MG_course/home/$USER/amplicon

cd $MY_PATH

for sample in `find /proj/g2017026/2017_MG_course/raw_data/amplicon/ -name "*_R1*.fastq"`;
do
  id=`basename ${sample%%_R1.fastq}`
  out_file=$MY_PATH/$id.fastq
  usearch10 -fastq_mergepairs $sample -fastqout $out_file --relabel "$id:"
  cat $out_file >> $MY_PATH/all_reads.fastq
done

### Question can you say something about the merging rate?! improve the rate maybe?!

usearch10 -fastq_filter all_reads.fastq -fastq_maxee 1 -fastaout all_reads.clean.fasta
usearch10 -fastx_uniques all_reads.clean.fasta -relabel Uniques_ -sizeout --minuniquesize 2 -fastaout all_reads.uniques.fasta
usearch10 -cluster_otus all_reads.uniques.fasta  -relabel OTUs_ -otus  all_reads.OTUs.fasta
usearch10 -usearch_global all_reads.fastq -db all_reads.OTUs.fasta -strand plus -id 0.97 -otutabout OTU_table.txt

### How do hit rates look for specific libs?!

DB=/proj/g2017026/2017_MG_course/data/usearch/rdp_16s_v16_sp.udb

usearch10 -sintax all_reads.OTUs.fasta -db $DB -tabbedout taxonomy.tax -strand both -sintax_cutoff 0.8
usearch10 -unoise3 all_reads.uniques.fasta -zotus all_reads.zOTUs.fasta
usearch10 -usearch_global all_reads.fastq -db all_reads.zOTUs.fasta -strand plus -id 0.97 -otutabout zOTU_table.txt
usearch10 -sintax all_reads.zOTUs.fasta -db ~/Data/tax_dbs/rdp_16s_v16_sp.udb -tabbedout ztaxonomy.tax -strand both -sintax_cutoff 0.8

Rscript ~/MG_course/2017_MG_course/scripts/ampliconplot.R
