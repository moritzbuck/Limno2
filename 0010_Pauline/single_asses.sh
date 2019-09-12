
mkdir $base_path/2000_assemblies/2100_single_samples

cd $base_path/2000_assemblies/2100_single_samples
mkdir $base_path/2000_assemblies/2100_single_samples/2110_megahit
cd $base_path/2000_assemblies/2100_single_samples/2110_megahit

#### Running sample-wise megahit assemblies

for s in ${samples[@]}
do
echo "Assembling $s with megahit"
megahit -1 $base_path/1000_processed_reads/1200_cleaned_reads/${s}_1P.gz -2 $base_path/1000_processed_reads/1200_cleaned_reads/${s}_2P.gz \
-r $base_path/1000_processed_reads/1200_cleaned_reads/${s}_1U.gz,$base_path/1000_processed_reads/1200_cleaned_reads/${s}_2U.gz \
-t 11 -o $s --out-prefix megahit_$s --continue
rm -r megahit_$s/intermediate_contigs
cp $s/megahit_$s.contigs.fa $base_path/2000_assemblies/2001_all_assemblies/
done



mkdir $base_path/2000_assemblies/2100_single_samples/2120_normalised_megahit
cd $base_path/2000_assemblies/2100_single_samples/2120_normalised_megahit

#### Running sample-wise megahit assemblies with the normalized reads

for s in ${samples[@]}
do
echo "Assembling $s with megahit"
megahit --12 $base_path/1000_processed_reads/1200_cleaned_reads/${s}_normed_paired.fastq.gz \
-t 11 -o $s --out-prefix normed_megahit_$s --continue
rm -r normed_megahit_$s/intermediate_contigs
cp $s/normed_megahit_$s.contigs.fa $base_path/2000_assemblies/2001_all_assemblies/
done

mkdir $base_path/2000_assemblies/2100_single_samples/2130_metaspades/
cd $base_path/2000_assemblies/2100_single_samples/2130_metaspades/

#### Running sample-wise metaspades assemblies

for s in ${samples[@]}
do

  echo "Submitting metaspades assembly $s"

  exec='#!/bin/bash'"\n

  module load spades/3.10.1\n
\n
  cp $base_path/1000_processed_reads/1200_cleaned_reads/${s}_clean.fastq /scratch/reads.fastq\n
  metaspades.py -m 500 --12 /scratch/reads.fastq -o /scratch/$s -t 15\n
  cp -r /scratch/$s $base_path/2000_assemblies/2100_single_samples/2130_metaspades/\n
  cp $s/scaffolds.fasta $base_path/2000_assemblies/2001_all_assemblies/spades_$s.contigs.fa\n
  "

  echo -e $exec | default_sbatch -o $base_path/2000_assemblies/2100_single_samples/2130_metaspades/$s.out -e $base_path/2000_assemblies/2100_single_samples/2130_metaspades/$s.err

done

cd $base_path/2000_assemblies/2100_single_samples/2130_metaspades/
for s in ${samples[@]}
do
  cp $s/scaffolds.fasta $base_path/2000_assemblies/2001_all_assemblies/metaspades_$s.contigs.fa
done


mkdir $base_path/2000_assemblies/2100_single_samples/2140_normed_metaspades/
cd $base_path/2000_assemblies/2100_single_samples/2140_normed_metaspades/

#### Running sample-wise normed metaspades assemblies

for s in ${samples[@]}
do

  echo "Submitting normed metaspades assembly $s"

  exec='#!/bin/bash'"\n

  module load spades/3.10.1\n
\n
  cp $base_path/1000_processed_reads/1200_cleaned_reads/${s}_normed_paired.fastq.gz /scratch/reads.fastq.gz\n
  unpigz /scratch/reads.fastq.gz \n
  metaspades.py -m 500 --12 /scratch/reads.fastq -o /scratch/$s -t 16\n
  cp -r /scratch/$s $base_path/2000_assemblies/2100_single_samples/2140_normed_metaspades/\n
  cp $s/scaffolds.fasta $base_path/2000_assemblies/2001_all_assemblies/normed_spades_$s.contigs.fa\n
  "

  echo -e $exec | default_sbatch -o $base_path/2000_assemblies/2100_single_samples/2140_normed_metaspades/$s.out -e $base_path/2000_assemblies/2100_single_samples/2140_normed_metaspades/$s.err -J normed_spades_$s

done
