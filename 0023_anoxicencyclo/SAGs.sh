#### SETTING UP env variables

base_path=/home/moritz/people/0023_anoxicencyclo/9997_MAGs
trimmomatic_path=/sw/apps/bioinfo/trimmomatic/0.36/rackham/
raws_folder=$base_path/0000_raws
python_hjelper=$base_path/helpers.py
phylophlan_path=$HOME/repos/github/phylophlan/

cd  $base_path

mkdir 1000_fastqc
mkdir 1000_fastqc/1100_prefiltering

#### Running fastQC

fastqc -o 1000_fastqc/1100_prefiltering $raws_folder/*/*.fastq.gz


#### Running trimmomatic

for f in `ls $raws_folder/*/*_R1_001.fastq.gz`
do
  echo "Processing $f"
  r=`echo $f | sed 's/_R1_/_R2_/'`
  out=`echo $f | sed 's/_L001_R1_001.fastq.gz//'`
  java -jar $trimmomatic_path/trimmomatic-0.36.jar PE -phred33 $f $r -baseout $out ILLUMINACLIP:$trimmomatic_path/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


#### moving the cleaned reads and cleaning up the names a bit

mkdir 2000_cleaned_reads
mv $raws_folder/*/*_[12][UP] 2000_cleaned_reads

samples=($(ls $base_path/2000_cleaned_reads | grep _1P | cut -d"_" -f1,2 | sort | uniq))

cd $base_path/2000_cleaned_reads
for s in ${samples[@]}
do
mv ${s}_1P ${s}_1P.fastq
mv ${s}_2P ${s}_2P.fastq
cat ${s}_1U ${s}_2U > ${s}_U.fastq
done
rm *[12]U

cd $base_path
mkdir 3000_assemblies
cd 3000_assemblies
mkdir $base_path/3000_assemblies/3001_all_assemblies/
#### Running sample-wise meta-spades

for s in ${samples[@]}
do
echo "Assembling $s with SC-spades"
  spades.py  --sc -1 $base_path/2000_cleaned_reads/${s}_1P.fastq -2 \
  $base_path/2000_cleaned_reads/${s}_2P.fastq -s \
  $base_path/2000_cleaned_reads/${s}_U.fastq -o $s -t 20 -k21,33,55,77,99,111 > /dev/null
  cp $s/scaffolds.fasta $base_path/3000_assemblies/3001_all_assemblies/$s.fasta
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
mkdir 4000_prokka

cd 4000_prokka
for f in ${samples[@]}${samples[@]}
do
  echo "Prokka on $f"
  prokka --outdir $f --force --prefix $f  --cpus 20 ../3000_assemblies/$f.10000bp.fasta  2> /dev/null
  sed -i 's/*//g' $f/$f.faa
  cp $f/$f.faa $f/SAG_$f.faa
done

cd $phylophlan_path
mkdir input/0023_SAGs

cp $base_path/4000_prokka/*/*.faa input/0023_SAGs
cp input/default_genomes/* input/0023_SAGs
python2.7 phylophlan.py --nproc 20  -i -t  0023_SAGs

cd $base_path
mkdir 5000_phylophlan
cd $base_path/5000_phylophlan
cp $phylophlan_path/output/0023_SAGs/* .

python $python_hjelper rename_tree 0023_SAGs.tree.int.nwk 0023_SAGs.named.tree -1

mkdir $base_path/0001_deliverables
cd $base_path/0001_deliverables
for f in ${samples[@]}
do
  mkdir SAG_$f
  cp ../4000_prokka/$f/$f.faa SAG_$f/proteins_$f.fasta
  cp ../4000_prokka/$f/$f.fna SAG_$f/genome_over10kb_$f.fasta
  cp ../4000_prokka/$f/$f.gff SAG_$f/annotation_$f.gff
  cp ../3000_assemblies/$f.pdf SAG_$f/assembly_plot_$f.pdf
  cp ../3000_assemblies/3001_all_assemblies/$f.fasta SAG_$f/all_scaffolds_$f.fasta
done

cp ../3000_assemblies/checkm.txt all_completnesses.txt
cat ../3000_assemblies/checkm.txt  | grep -v 1000bp | grep -v 2500bp | sed 's/.10000bp//' > only_10kb_completnesses.txt
cp ../5000_phylophlan/0023_SAGs.named.tree  phylophlan.tree
cp ~/Dropbox/Documents/Uppsala/Limno/0001_Chloroflexi_SAGs/script.sh this_script.sh
cp $python_hjelper python_hjelper_script.py
cp ~/Dropbox/Documents/Uppsala/Limno/0001_Chloroflexi_SAGs/README.txt .

cd ..
tar czvf SAGs_and_info.tar.gz 0001_deliverables/



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
