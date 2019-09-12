#!/bin/base

#### SETTING UP env variables


base_path=/home/moritz/bils_folder/people/0010_Pauline
if [ $HOSTNAME = "omneya-desktop" ]
then
  base_path=/home/moritz/people/0010_Pauline
fi


trimmomatic_path=/home/moritz/share/Trimmomatic-0.36/
raws_folder=$base_path/0000_raws/0100_reads
python_hjelper=$HOME/Dropbox/Documents/Uppsala/Limno/0010_Pauline/helpers.py


alias default_sbatch="sbatch -A b2011035 -p node -n 16 -t 10-00:00:00 --mail-type=ALL -C mem512GB --mail-user murumbii@gmail.com"

samples=($(find $raws_folder/* -maxdepth 0  -type d  | grep -v Report| cut -d'/' -f9))
if [ $HOSTNAME = "omneya-desktop" ]
then
  samples=($(find $raws_folder/* -maxdepth 0  -type d  | grep -v Report| cut -d'/' -f8))
fi

num_threads=10

cd  $base_path

#### Running fastQC

mkdir 1000_processed_reads
mkdir 1000_processed_reads/1100_fastqc

fastqc -o 1000_processed_reads/1100_fastqc `find $raws_folder -name "*.fastq.gz"`

#### Running trimmomatic
mkdir 1000_processed_reads/1200_cleaned_reads

for s in ${samples[@]}
do
  echo "Processing $s"
  f=`find $raws_folder/$s/02-FASTQ/ -name  "*_R1_*"`
  r=`find $raws_folder/$s/02-FASTQ/ -name  "*_R2_*"`

  out=$base_path/1000_processed_reads/1200_cleaned_reads/$s

  java -jar $trimmomatic_path/trimmomatic-0.36.jar PE -phred33 $f $r -baseout $out ILLUMINACLIP:$trimmomatic_path/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

cd $base_path/1000_processed_reads/1200_cleaned_reads/
#### noirmalizing the reads with the digital normalization function from the khmer software package

for s in ${samples[@]}
do
  pigz ${s}_1P  ${s}_2P  ${s}_2U ${s}_1U
  interleave-reads.py ${s}_1P.gz ${s}_2P.gz -o ${s}_clean.fastq.gz
  normalize-by-median.py -p -k 20 -C 20 -N 4 -x 8e9 -s $s.kh -o ${s}_norm_step_1.fastq ${s}_clean.fastq
  filter-abund.py --threads 11 -V ${s}.kh ${s}_norm_step_1.fastq -o ${s}_norm_step_2.fastq
  normalize-by-median.py -C 5 -k 20 -N 4 -x 8e9  -o ${s}_normed.fastq ${s}_norm_step_2.fastq
  extract-paired-reads.py -p ${s}_normed_paired.fastq -s /dev/null ${s}_normed.fastq
  rm ${s}_norm_step_1.fastq ${s}_norm_step_2.fastq ${s}.kh  ${s}_normed_paired.fastq
  pigz ${s}_normed_paired.fastq  ${s}_clean.fastq
done

cd $base_path
#### computiong mash distances for reads
for s in ${samples[@]}
do
  mash sketch -r -p 11 -k 21 -s 10000  $base_path/1000_processed_reads/1200_cleaned_reads/${s}_clean.fastq.gz
done

mkdir $base_path/1000_processed_reads/1300_mash
cd $base_path/1000_processed_reads/1300_mash
mv $base_path/1000_processed_reads/1200_cleaned_reads/*.msh .

echo > distances.tsv
ls *.msh| sed 's/.msh//' >> distances.tsv

for f in `ls *.msh`
do
  ass=${f%%.msh}
  echo $ass > temp.tsv
  mash dist $f *.msh | cut -f3 >> temp.tsv
  paste distances.tsv temp.tsv > temp2.tsv
  mv temp2.tsv distances.tsv
done
rm temp.tsv

python $python_hjelper mash_dist_plot

mkdir $base_path/1000_processed_reads/1400_taxonomy
mkdir $base_path/1000_processed_reads/1400_taxonomy/1410_kaiju/
mkdir $base_path/1000_processed_reads/1400_taxonomy/1410_kaiju/1411_proGenome
cd $base_path/1000_processed_reads/1400_taxonomy/1410_kaiju/1411_proGenome


db_path=~/DataBases/kaiju/proGenome/

for s in ${samples[@]}
do
  echo "Taxonomically annotating $s with kaiju"
  cd $base_path/1000_processed_reads/1400_taxonomy/1410_kaiju/1411_proGenome

  if [ ! -d $s ]
  then
    mkdir $s
  fi

  if [ ! -f ${s}/kaiju.out ]
  then
    kaiju -t $db_path/nodes.dmp -f $db_path/kaiju_db.fmi -i $base_path/1000_processed_reads/1200_cleaned_reads/${s}_1P.gz -j \
    $base_path/1000_processed_reads/1200_cleaned_reads/${s}_2P.gz -o ${s}/kaiju.out -z $num_threads
    kaiju2krona -u -t $db_path/nodes.dmp -n $db_path/names.dmp -i ${s}/kaiju.out -o ${s}/kaiju.out.krona
    ktImportText -o ${s}/$s.html ${s}/kaiju.out.krona
    kaijuReport -t $db_path/nodes.dmp -n $db_path/names.dmp -i $s/kaiju.out -r family -o $s/kaiju.out.summary
  fi

done

db_path=~/DataBases/kaiju/nr/
mkdir $base_path/1000_processed_reads/1400_taxonomy/1410_kaiju/1412_nr
for s in ${samples[@]}
do
  echo "Taxonomically annotating $s with kaiju"
  cd $base_path/1000_processed_reads/1400_taxonomy/1410_kaiju/1412_nr

  if [ ! -d $s ]
  then
    mkdir $s
  fi

  if [ ! -f ${s}/kaiju.out ]
  then
    kaiju -t $db_path/nodes.dmp -f $db_path/kaiju_db_nr_euk.fmi -i $base_path/1000_processed_reads/1200_cleaned_reads/${s}_clean.fastq -o ${s}/kaiju.out -z $num_threads
    kaiju2krona -u -t $db_path/nodes.dmp -n $db_path/names.dmp -i ${s}/kaiju.out -o ${s}/kaiju.out.krona
    ktImportText -o ${s}/$s.html ${s}/kaiju.out.krona
    kaijuReport -t $db_path/nodes.dmp -n $db_path/names.dmp -i $s/kaiju.out -r family -o $s/kaiju.out.summary
  fi
done






cd $base_path

mkdir 2000_assemblies
cd 2000_assemblies
mkdir 2001_all_assemblies

./single_asses.sh

./coasses.sh

#### Filtering the assemblies using the custom python script

cd $base_path/2000_assemblies/2001_all_assemblies/
for f in `ls *.fa | grep -v ".le"`
do
  echo "Filtering $f"
  out_f=`echo $f | sed 's/.contigs./.le2500bp.contigs./'`
  if [ ! -f $out_f  ]; then
  python $python_hjelper filter_fasta_by_len 2500 $f $out_f
  fi
  out_f=`echo $f | sed 's/.contigs./.le1000bp.contigs./'`
  if [ ! -f $out_f  ]; then
    python $python_hjelper filter_fasta_by_len 1000 $f $out_f
  fi
done


cd $base_path
mkdir 4000_maps

#### Mapping subsets of all libraries to all assemblies

echo "making bbmap mapping rates"
mkdir $base_path/4000_maps/4100_bbmap
cd $base_path/4000_maps/4100_bbmap
for f in `ls $base_path/2000_assemblies/2001_all_assemblies/*.fa`
do

  cd $base_path/4000_maps/4100_bbmap
  ass=`basename $f`
  ass=${ass%%.contigs.fa}
  echo "Doing $ass"
  if [ ! -d $ass  ]; then
    mkdir $ass
  fi

  cd $ass

  if [ ! -d ref  ]; then
    bbmap.sh ref=$f 2> /dev/null
  fi

  for s in ${samples[@]}
  do
    if [ ! -f $s.rates  ]; then
      echo "Mapping $s"
      bbmap.sh in=$base_path/1000_processed_reads/1200_cleaned_reads/${s}_1P.gz in2=$base_path/1000_processed_reads/1200_cleaned_reads/${s}_2P.gz reads=200000 threads=10 2> $s.rates
    fi
  done
done


mkdir $base_path/500_assemblystats
cd $base_path/500_assemblystats

## extract mapping rates from the bbmap outputs
grep mated $base_path/400_maps/*/*.rates | cut -f1,2 | sed 's/:mated pairs://' | sed 's/%//'  | sed 's/\//\t/' | sed 's/.rates//'  > all_mapping_rates.tsv

#### Running assemstats (from the khmer software package) on all assemblies to get basic assembly stats
assemstats.py 0 $base_path/2000_assemblies/2001_all_assemblies/*.fa > all_assembly_stats.txt

#### Using the python script to compile all data in a nice table
python $python_hjelper format_assembly_tables all_assembly_stats.txt all_mapping_rates.tsv all_assemblies_table.csv

#### copying main files to a deliverable folder

mkdir $base_path/999_deliverable

cp $HOME/Dropbox/Documents/Uppsala/NBIS/3726_mican/script.sh $base_path/999_deliverable/this_script.sh
cp $HOME/Dropbox/Documents/Uppsala/NBIS/3726_mican/helpers.py $base_path/999_deliverable/
cp $HOME/Dropbox/Documents/Uppsala/NBIS/3726_mican/3726* $base_path/999_deliverable/
