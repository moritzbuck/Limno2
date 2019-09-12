#### SETTING UP env variables


source ~/people/0010_Pauline/arctic/bin/activate

base_path=$HOME/people/0006_Erken_Genome_Project
raws_folder=$HOME/people/0006_Erken_Genome_Project/0000_rawdata/0100_SAG_data/0110_libtest/
python_hjelper=$HOME/temp/helpers.py
phylophlan_path=$HOME/repos/github/phylophlan/
THREADS=20
clean_reads=1000_reads_proc/1200_trimmomatic/
nb_reads_HiX=437534089
declare -a rare_levs=(500 250 100 50 0)

cd  $base_path

samples=($(ls $clean_reads  | cut -f1 -d"_" | sort | uniq))
models=`ls $clean_reads  | cut -f1 -d"_" | sort | uniq | grep QFX | cut -f 1-5 -d"-"`


temp_path=1000_reads_proc/1100_fastqc
mkdir -p $temp_path/1110_beforefilter

#### Running fastQC
module load bioinfo-tools
module load FastQC
module load seqtk
module load trimmomatic
module load bbmap

fastqc --threads 6 -o $temp_path/1110_beforefilter $raws_folder/*/*.fastq.gz

#### Running trimmomatic


for f in `ls $raws_folder/*-NXT/*_R1_001.fastq.gz`
do
  echo "Processing $f"
  r=`echo $f | sed 's/_R1_/_R2_/'`
  out=`echo $f | sed 's/_L001_R1_001//'`
  cp $f /scratch/
  cp $r /scratch/
  java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -phred33 /scratch/`basename $f` /scratch/`basename $r` -baseout /scratch/`basename $out` -threads $THREADS ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> /scratch/`basename $out`.log
  cp /scratch/`basename ${out%%.fastq.gz}`*  `dirname $out`
done

for f in `ls $raws_folder/*/*_R1_001.fastq.gz | grep -v NXT`
do
  r=`echo $f | sed 's/_R1_/_R2_/'`
  out=`echo $f | sed 's/_L[0-9]*_R1_001//'`
  fin=`basename ${out%%.fastq.gz}`_1P.fastq.gz

  if [ ! -f $clean_reads/$fin  ]
  then
      echo "Processing $f"

      cp $f /scratch/
      cp $r /scratch/
      java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -phred33 /scratch/`basename $f` /scratch/`basename $r` -baseout /scratch/`basename $out` -threads $THREADS ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> /scratch/`basename $out`.log
      cp /scratch/`basename ${out%%.fastq.gz}`*  $clean_reads
  fi
done

mkdir -p $clean_reads
mv $raws_folder/*/*_[12][UP].fastq.gz $clean_reads





#### moving the cleaned reads and cleaning up the names a bit

samples=($(ls $clean_reads  | cut -f1 -d"_" | sort | uniq))
models=`ls $clean_reads  | cut -f1 -d"_" | sort | uniq | grep QFX | cut -f 1-5 -d"-"`

for nb_cells in ${rare_levs[@]}
do
  for s in ${samples[@]}
  do
    run=0
    #### Running sample-wise meta-spades
    if [ $nb_cells = 0 ]
    then
      ass_name="raw"
      run=1
    else
      ass_name=${nb_cells}cells
      run=$(echo $models | grep -c  `echo $s | cut -f 1-5 -d"-"`)
    fi
    if [ $run = 1 ]
    then
      mkdir -p 3000_assemblies/$s/
      ass_path=3000_assemblies/$s/3010_$ass_name
      fwd=`ls $clean_reads | grep $s | grep _1P`
      rev=`ls $clean_reads | grep $s | grep _2P`
      if [ ! -f $ass_path/scaffolds.fasta ]
    	then
        if [ $nb_cells = 0 ]
        then
           zcat $clean_reads/$fwd >  /scratch/${fwd%%.gz}
           zcat $clean_reads/$rev >  /scratch/${rev%%.gz}
        else
          seqtk sample -s100 $clean_reads/$fwd `echo $nb_reads_HiX/$nb_cells | bc ` >  /scratch/${fwd%%.gz}
          seqtk sample -s100 $clean_reads/$rev `echo $nb_reads_HiX/$nb_cells | bc ` >  /scratch/${rev%%.gz}
        fi
        echo "Assembling $s as $ass_name with SPAdes"
        spades.py --sc -1 /scratch/${fwd%%.gz} -2 /scratch/${rev%%.gz} -o $ass_path -t 20 > 3000_assemblies/$s/$ass_name.log
      fi
      if [ ! -f $ass_path/$s.pdf ]
      then
        echo "QC assembly for $s"
        python $python_hjelper spades_assembly_qc_plot $ass_path/scaffolds.fasta $ass_path/$s.pdf 2> /dev/null
      fi
    fi
  done
done

#### Filtering
for nb_cells in ${rare_levs[@]}
do
  for s in ${samples[@]}
  do
    run=0
    #### Running sample-wise meta-spades
    if [ $nb_cells = 0 ]
    then
      ass_name="raw"
      run=1
    else
      ass_name=${nb_cells}cells
      run=$(echo $models | grep -c  `echo $s | cut -f 1-5 -d"-"`)
    fi
    if [ $run = 1 ]
    then
      ass_path=3000_assemblies/$s/3010_$ass_name
      out_path=3000_assemblies/$s/3011_filtered/
      if  [  -f $ass_path/scaffolds.fasta ]
      then
        lengths=(1000 2500 10000)
        for len in ${lengths[@]}
        do
          echo "Filtering $s with $ass_name at ${len}bp"

          mkdir -p $out_path
          out_p=$out_path/${s}_${ass_name}_${len}bp.fasta
          python $python_hjelper filter_fasta_by_len $len $ass_path/scaffolds.fasta $out_p
        done
      fi
    fi
  done
done
mkdir -p 4000_assembly_analysis/4100_checkm/
mkdir -p 4000_assembly_analysis/4999_genomes/

cp `find 3000_assemblies -name "*bp.fasta" `  4000_assembly_analysis/4999_genomes/

checkm lineage_wf -t 10 -x fasta  4000_assembly_analysis/4999_genomes/ 4000_assembly_analysis/4100_checkm/data > 4000_assembly_analysis/4100_checkm/checkm.txt

mkdir -p 3000_assemblies/used_assemblies
cp `find 3000_assemblies/ -name "*raw_*bp.fasta"  | grep -f <(echo $models  | tr ' ' $'\n') `  3000_assemblies/used_assemblies/




####################### Mapping for other stuff ############
model_list=($(echo $models | grep -v NXT))
module load bbmap

root=`pwd`

for mod in ${model_list[@]}
do
  echo doing $mod

#  libs=($(ls $clean_reads  | cut -f1 -d"_" | sort | uniq | grep $mod | grep -v NXT))

  cat /dev/null > /scratch/in_reads_fwd.fastq
  cat /dev/null > /scratch/in_reads_rev.fastq
  cat /dev/null > /scratch/out_reads_fwd.fastq
  cat /dev/null > /scratch/out_reads_rev.fastq

  echo "==================  out fwd =============="
  for ll in `find $root/1000_reads_proc/1200_trimmomatic/ | grep -v $mod | grep _1P | grep -v NXT`
  do
    unpigz -c -p 20 $ll | head -n 1000000 >>  /scratch/out_reads_fwd.fastq
  done
  echo "==================  out rev =============="
  for ll in `find $root/1000_reads_proc/1200_trimmomatic/ | grep -v $mod | grep _2P | grep -v NXT`
  do
    unpigz -c -p 20 $ll | head -n 1000000 >>  /scratch/out_reads_rev.fastq
  done

  for ref in `find 3000_assemblies/used_assemblies/ | grep $mod | grep -v NXT`
  do
    ff=`basename $ref`
    pat=4000_assembly_analysis/4300_mapping/$ff
    mkdir -p $pat
    cd $pat
    cell=`echo $ff | cut -d'-' -f1-3`
    if [ ! -f out_map.log ]
    then
        echo this ass $ref
        bbmap.sh ref=$root/$ref in=/scratch/out_reads_fwd2.fastq in2=/scratch/out_reads_rev2.fastq 2> out_map.log
    fi
    cd $root
  done
done


model_list=($(echo $models))

for mod in ${model_list[@]}
do
  echo doing $mod

#  libs=($(ls $clean_reads  | cut -f1 -d"_" | sort | uniq | grep $mod | grep -v NXT))

  cat /dev/null > /scratch/out_reads_fwd2.fastq
  cat /dev/null > /scratch/out_reads_rev2.fastq
  cat /dev/null > /scratch/in_reads_fwd2.fastq
  cat /dev/null > /scratch/in_reads_rev2.fastq

  echo "==================  out fwd =============="
  for ll in `find $root/1000_reads_proc/1200_trimmomatic/ | grep -v $mod | grep _1P | grep NXT`
  do
    unpigz -c -p 20 $ll | head -n 1000000 >>  /scratch/out_reads_fwd2.fastq
  done
  echo "==================  out rev =============="
  for ll in `find $root/1000_reads_proc/1200_trimmomatic/ | grep -v $mod | grep _2P | grep NXT`
  do
    unpigz -c -p 20 $ll | head -n 1000000 >>  /scratch/out_reads_rev2.fastq
  done

  echo "==================  out fwd =============="
  for ll in `find $root/1000_reads_proc/1200_trimmomatic/ | grep  $mod | grep _1P | grep NXT`
  do
    unpigz -c -p 20 $ll | head -n 6000000 >>  /scratch/in_reads_fwd2.fastq
  done
  echo "==================  out rev =============="
  for ll in `find $root/1000_reads_proc/1200_trimmomatic/ | grep  $mod | grep _2P | grep NXT`
  do
    unpigz -c -p 20 $ll | head -n 6000000 >>  /scratch/in_reads_rev2.fastq
  done


  for ref in `find 3000_assemblies/used_assemblies/ | grep $mod | grep NXT`
  do
    ff=`basename $ref`
    pat=4000_assembly_analysis/4300_mapping/$ff
    mkdir -p $pat
    cd $pat
    cell=`echo $ff | cut -d'-' -f1-3`
    if [ ! -f out_map.log ]
      then
        echo this ass $ref
        bbmap.sh ref=$root/$ref in=/scratch/out_reads_fwd2.fastq in2=/scratch/out_reads_rev2.fastq 2> out_map.log
      fi
      cd $root
  done
done



grep mapped 4000_assembly_analysis/4300_mapping/*/out_map.log | sed 's#.*mapping/\(RE.*\).fasta.*mapped:\(.*\)#\1\t\2#' | awk '{print $1 "\t" $2}' | tr -d '%' > 4000_assembly_analysis/4300_mapping/mapped.txt
grep mated 4000_assembly_analysis/4300_mapping/*/out_map.log | sed 's#.*mapping/\(RE.*\).fasta.*mated pairs:\(.*\)#\1\t\2#' | awk '{print $1 "\t" $2}' | tr -d '%'  > 4000_assembly_analysis/4300_mapping/mated.txt
