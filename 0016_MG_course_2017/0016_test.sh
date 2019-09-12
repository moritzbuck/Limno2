module load bioinfo-tools
module load FastQC
module load trimmomatic
module load megahit




sample=feces_004
base=/home/moritz/MG_course/2017_MG_course/
kaiju_db=$base/data/kaijudb/proGenomes/
opat=$base/../home/moritz/shotgun/
export PATH=$PATH:$base/bin::~/.local/bin/:~/.local/bin

mkdir $opat/trimmomatic/$sample

fastqc -o fastqc /home/moritz/MG_course/raw_data/shotgun/feces_004_R1.fastq

java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -phred33 -threads 8  $base/raw_data/shotgun/${sample}_R1.fastq $base/raw_data/shotgun/${sample}_R2.fastq  -baseout $opat/trimmomatic/${sample}.fastq ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# up to 8 min

mash sketch -p 20 -r $opat/trimmomatic/${sample}/${sample}_1P.fastq -o $opat/mash/${sample}
megahit -1 $opat/trimmomatic/${sample}/${sample}_1P.fastq -2 $opat/trimmomatic/${sample}/${sample}_2P.fastq -o $opat/megahit/${sample}
mkdir $opat/kaiju/${sample}
kaiju -t $kaiju_db/nodes.dmp -f ${kaiju_db}/kaijudb.fmi -i $opat/trimmomatic/${sample}/${sample}_1P.fastq -j $opat/trimmomatic/${sample}/${sample}_2P.fastq -o $opat/kaiju/${sample}/$sample -z 10
#up to 12ish minutes

mkdir $opat/kaiju/${sample}
mkdir $opat/humann/${sample}

humann2 --input ../subset/feces_004_sub_1P.fastq

bbmap.sh ref=ASS.fasta in=FWD_LIB in2=REV_LIB out=LIB_2_ASS.sam bamscript=LIB_2_ASS.sam2ban.sh 
