#### SETTING UP env variables

base_path=/home/moritz/Data/people/0001_Chloroflexi_SAGs
trimmomatic_path=/home/moritz/share/Trimmomatic-0.36/
raws_folder=/home/moritz/Data/raws/0001_Chloroflexi_SAGs
python_hjelper=$HOME/Dropbox/Documents/Uppsala/NBIS/3726_mican/helpers.py
phylophlan_path=$HOME/repos/mercurial/phylophlan/

cd  $base_path

mkdir 100_fastqc
mkdir 100_fastqc/110_prefiltering

#### Running fastQC

fastqc -o 100_fastqc/110_prefiltering $raws_folder/*/*.fastq.gz


#### Running trimmomatic

for f in `ls $raws_folder/*/*_R1_001.fastq.gz`
do
  echo "Processing $f"
  r=`echo $f | sed 's/_R1_/_R2_/'`
  out=`echo $f | sed 's/_L001_R1_001.fastq.gz//'`
  java -jar $trimmomatic_path/trimmomatic-0.36.jar PE -phred33 $f $r -baseout $out ILLUMINACLIP:$trimmomatic_path/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

mkdir 200_cleaned_reads
mv $raws_folder/*/*_[12][UP] cd/


#### moving the cleaned reads and cleaning up the names a bit

samples=($(ls $base_path/200_cleaned_reads | grep _1P | cut -d"_" -f1,2 | sort | uniq))

cd $base_path/200_cleaned_reads
for s in ${samples[@]}
do
mv ${s}_1P ${s}_1P.fastq
mv ${s}_2P ${s}_2P.fastq
cat ${s}_1U ${s}_2U > ${s}_U.fastq
done

cd $base_path
mkdir 300_assemblies
cd 300_assemblies

#### Running sample-wise meta-spades

for s in ${samples[@]}
do
echo "Assembling $s with metaspades"
  spades.py  --sc -1 $base_path/200_cleaned_reads/${s}_1P.fastq -2 \
  $base_path/200_cleaned_reads/${s}_2P.fastq -s \
  $base_path/200_cleaned_reads/${s}_U.fastq -o $s -t 1 > /dev/null
  cp $s/megahit_$s.contigs.fa $base_path/300_assemblies/301_all_assemblies/
done

#### Custom spades assembly QC plots (using the python script)

for s in ${samples[@]}${samples[@]}
do
  echo "QC assembly for $s"
  python $python_hjelper spades_assembly_qc_plot $s/scaffolds.fasta $s.pdf 2> /dev/null
done

#### Filtering

for f in ${samples[@]}${samples[@]}
do
  echo "Filtering $f"
  out_f=$f.1000bp.fasta
  python $python_hjelper filter_fasta_by_len 1000 $f/scaffolds.fasta $out_f
  out_f=$f.2500bp.fasta
  python $python_hjelper filter_fasta_by_len 2500 $f/scaffolds.fasta $out_f
  out_f=$f.10000bp.fasta
  python $python_hjelper filter_fasta_by_len 10000 $f/scaffolds.fasta $out_f
done

checkm lineage_wf -t 10 -x fasta . checkm > checkm.txt

cd $base_path
mkdir 400_prokka

cd 400_prokka
for f in ${samples[@]}${samples[@]}
do
  echo "Prokka on $f"
  prokka --outdir $f --force --prefix $f  --cpus 10 ../300_assemblies/$f.10000bp.fasta  2> /dev/null
  sed -i 's/*//g' $f/$f.faa
  cp $f/$f.faa $f/SAG_$f.faa
done

cd $phylophlan_path
mkdir input/0001_Chloroflexi_SAGs

cp $base_path/400_prokka/*/SAG_*.faa input/0001_Chloroflexi_SAGs
cp input/default_genomes/* input/0001_Chloroflexi_SAGs
python phylophlan.py --nproc 16  -i -t  0001_Chloroflexi_SAGs

cd $base_path
mkdir 500_phylophlan
cd $base_path/500_phylophlan
cp $phylophlan_path/output/0001_Chloroflexi_SAGs/* .

python $python_hjelper rename_tree 0001_Chloroflexi_SAGs.tree.int.nwk 0001_Chloroflexi_SAGs.named.tree -1

mkdir $base_path/001_deliverables
cd $base_path/001_deliverables
for f in ${samples[@]}${samples[@]}
do
  mkdir SAG_$f
  cp ../400_prokka/$f/$f.faa SAG_$f/proteins_$f.fasta
  cp ../400_prokka/$f/$f.fna SAG_$f/genome_over10kb_$f.fasta
  cp ../400_prokka/$f/$f.gff SAG_$f/annotation_$f.gff
  cp ../300_assemblies/$f.pdf SAG_$f/assembly_plot_$f.pdf
  cp ../300_assemblies/$f/scaffolds.fasta SAG_$f/all_scaffolds_$f.fasta
done

cp ../300_assemblies/checkm.txt all_completnesses.txt
cat ../300_assemblies/checkm.txt  | grep -v 1000bp | grep -v 2500bp | sed 's/.10000bp//' > only_10kb_completnesses.txt
cp ../500_phylophlan/0001_Chloroflexi_SAGs.named.tree  phylophlan.tree
cp ~/Dropbox/Documents/Uppsala/Limno/0001_Chloroflexi_SAGs/script.sh this_script.sh
cp $python_hjelper python_hjelper_script.py
cp ~/Dropbox/Documents/Uppsala/Limno/0001_Chloroflexi_SAGs/README.txt .

cd ..
tar czvf SAGs_and_info.tar.gz 001_deliverables/



####################### Mapping for other stuff ############

cd $base_path/400_maps
for f in `ls $base_path/300_assemblies/301_all_assemblies/spades*.le2500bp*.fa`
do
  cd $base_path/400_maps
  ass=`basename $f`
  ass=${ass%%.contigs.fa}
  echo "Doing $ass"
  mkdir $ass
  cd $ass
  for s in ${samples[@]}
  do
    echo "Full Mapping $s"
    bbmap.sh  in=$base_path/200_cleaned_reads/${s}_L001_1P.gz showprogress=333333 in2=$base_path/200_cleaned_reads/${s}_L001_2P.gz out=$s.sam threads=10 2> $s.full.rates
    sambamba_v0.6.6  view  -t 10 -S -h -f bam  $s.sam -o $s.bam
    sambamba_v0.6.6  sort  -t 10 -o ${s}_sorted.bam $s.bam
    sambamba_v0.6.6 index ${s}_sorted.bam
  done
done
