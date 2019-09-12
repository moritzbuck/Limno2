
DIR=$HOME/people/0010_Pauline
cd
cd $DIR

QSCRIPT="sbatch -D $DIR -A snic2017-1-616 -t 5-00:00:00 -n 20 -p node --mail-user murumbii@gmail.com --mail-type=ALL"

all_samples=`find -H 0000_raws/0100_reads/genomic/*  -type "d" | grep -v FASTQ  | rev | cut -f1 -d"/" | rev | sed 's/^/1000_processed_reads\//' | sed 's/$/\/done/' `
snakemake  --cores 20 $all_samples
all_single_asses=` find -H 0000_raws/0100_reads/genomic/*  -type "d" | grep -v FASTQ  | rev | cut -f1 -d"/" | rev | sed 's/^/1000_processed_reads\//' | sed 's/$/\/assemblies\/megahit\/mapping\/map_table.tsv/'`

all_megahit=`ls 0000_raws/0200_coasses | grep -v Full | sed 's/^/1500_assemblies\//' | sed 's/.txt$/\/megahit\//' `
all_spades=`ls 0000_raws/0200_coasses | grep -v Full | sed 's/^/1500_assemblies\//' | sed 's/.txt$/\/spades\//' `
snakemake  --cluster "$QSCRIPT" --jobs 500 --cores 20 $all_single_asses
snakemake  --cluster "$QSCRIPT" --cores 500 --local-cores 20
