#!/bin/sh
#SBATCH -o /home/moritz/0016_shotgun.out
#SBATCH -e /home/moritz/0016_shotgun.err

module load trimmomatic
trimmomatic_path=/home/moritz/share/Trimmomatic-0.36/
trimmomatic_path=$TRIMMOMATIC_HOME

export PATH=$PATH:`pwd`/bin

module load bbmap
module load MetaBat
module load samtools

cd ~/share/
wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar
tar xvf mash-Linux64-v2.0.tar
ln $HOME/share/mash-Linux64-v2.0/mash $HOME/bin/

mkdir /proj/g2017026/0016_MG_course_2017/shotgun/
mkdir /proj/g2017026/0016_MG_course_2017/all_assemblies/

cd /proj/g2017026/0016_MG_course_2017/shotgun/

find . -name "*_R1*" > all_fwds.txt

for f in `cat all_fwds.txt`;
do
  r=`echo $f | sed 's/_R1.fastq/_R2.fastq/'`
  out=`echo $f | sed 's/_R1.fastq/.fastq/'`
  java -jar $trimmomatic_path/trimmomatic-0.36.jar PE -phred33 $f $r -baseout $out ILLUMINACLIP:$trimmomatic_path/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


sed 's/_R1.fastq/_1P.fastq/' all_fwds.txt > clean_fwds.txt
sed 's/_R1.fastq/_1P.fastq.msh/' all_fwds.txt > mashes.txt
cp clean_fwds.txt assembly_files/all.txt

for f in `cat clean_fwds.txt`;
do
  mash sketch -p 20 -r $f
done

ls common/mashes/*.msh > mashes.txt
for f in `cat mashes.txt`;
do
  mash dist -p 20 -l $f mashes.txt
done > mash_distances.txt



for f in `cat all_fwds.txt`;
do
  r=`echo $f | sed 's/_R1.fastq/_R2.fastq/'`
  out=`echo $f | sed 's/_R1.fastq/.fastq/'`

  ft=`echo $f | sed 's/_R1.fastq/_1P.fastq/'`
  rt=`echo $f | sed 's/_R1.fastq/_2P.fastq/'`
  fu=`echo $f | sed 's/_R1.fastq/_1U.fastq/'`
  ru=`echo $f | sed 's/_R1.fastq/_2U.fastq/'`
  mo=`echo $f | sed 's/_R1.fastq//'`
  echo starting $out
  if [ ! -f ${mo}/final.contigs.fa ]
  then
    rm -r $mo
    megahit -1 $ft -2 $rt -r $fu,$ru -o $mo 2> /dev/null > /dev/null
  fi
  cp ${mo}/final.contigs.fa /proj/g2017026/0016_MG_course_2017/all_assemblies/`basename $mo`.single.fasta
done

rm -r `find . -name "intermediate_contigs"`

for subset in `ls assembly_files`
do
  name=`basename ${subset%%.txt}`
  fwds=`grep -f $subset clean_fwds.txt | tr '\n' ','`
  revs=`grep -f $subset clean_fwds.txt | sed 's/_1P.fastq/_2P.fastq/' | tr '\n' ','`
  singles=`grep -f $subset clean_fwds.txt | sed 's/_1P.fastq/_2U.fastq/' | tr '\n' ','``grep -f $subset clean_fwds.txt | sed 's/_1P.fastq/_1U.fastq/' | tr '\n' ','`
  megahit -1 ${fwds::-1} -2 ${revs::-1} -r ${singles::-1} -o $name
  cp $name/final.contigs.fa ../all_assemblies/$name.fasta
done

rm -r `find . -name "intermediate_contigs"`

mkdir maps
for ass in `ls /proj/g2017026/0016_MG_course_2017/all_assemblies/*.fasta`
do
  echo mapping $ass
  name=`basename $ass`
  mkdir maps/$name
  if [ ! -d maps/$name/ref ]
  then
    bbmap.sh ref=$ass path=maps/$name
  fi

  for f in `cat all_fwds.txt`;
  do
    r=`echo $f | sed 's/_R1.fastq/_R2.fastq/'`
    out=`echo $f | sed 's/_R1.fastq/.fastq/'`
    ft=`echo $f | sed 's/_R1.fastq/_R1.fastq/'`
    rt=`echo $f | sed 's/_R1.fastq/_R2.fastq/'`
    mo=`echo $f | sed 's/_R1.fastq//'`
    mo=`basename $mo`
    echo mapping $mo
    if [ ! -f maps/$name/${mo}_sorted.bam ]
    then
      bbmap.sh path=maps/$name reads=1000000 in=$ft in2=$rt out=maps/$name/$mo.sam bamscript=maps/$name/$mo.sh #2> maps/$name/$mo.txt
      maps/$name/$mo.sh
      rm maps/$name/$mo.sam
    fi
  done
done

for ass in `ls /proj/g2017026/0016_MG_course_2017/all_assemblies/*.fasta`
do
  echo binning $ass
  name=`basename $ass`
  runMetaBat.sh $ass `find maps/$name -name "*.bam"`
  mkdir  maps/$name/binning/
  mv ${name}* maps/$name/binning/
done

mkdir ../all_bins
for fold in `find . -name "*metabat-bins"`
do
  name=`echo $fold | cut -f 3 -d'/' `
  for b in `ls $fold  | grep bin.*fa`
  do
    cp $fold/$b ../all_bins/${name%%.fasta}.$b
  done
done

cd ../all_bins
checkm lineage_wf -t 20 -x fa . checkm > checkm.txt
