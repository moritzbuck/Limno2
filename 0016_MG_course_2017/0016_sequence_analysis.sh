
module load biopython
module load cd-hit


cd /proj/g2017026/0016_MG_course_2017/metagenome

python ../0016_fasta_filter.py  1000 all_contigs.fa larger_1000.fa
prokka --metagenome larger_1000.fa --outdir annotations --prefix combined

cd-hit-est -i annotations/combined.ffn -o annotations/combine_clusted.ffn -c 0.95 -T 16

mkdir mapps
for f in `find /proj/g2017026/0016_MG_course_2017/shotgun/ -name "*_R1*"`;
do
  r=`echo $f | sed 's/_R1.fastq/_R2.fastq/'`
  lib_name=`basename ${f%%_R1.fastq}`
  echo doing bbmap
  bbmap.sh ref=annotations/combine_clusted.ffn in=$f in2=$r out=mapps/$lib_name.sam 2> mapps/$lib_name.out
  sambamba view -p -S -t 20 -f bam mapps/$lib_name.sam -o temp.bam

done
