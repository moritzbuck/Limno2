base_path=/home/moritz/Data/people/0003_Vidy
raws_folder=$base_path/000_data/020_reads/
python_hjelper=$HOME/Dropbox/Documents/Uppsala/Limno/0002_stump/helper.py
usearch=usearch9
usearch_path=$base_path/200_usearch/
nproc=11

cd  $base_path

mkdir $base_path/100_formated_data
mkdir $base_path/100_formated_data/deplexed
for f in `ls $base_path/000_data/010_jsons/*.json`
do
  python $python_hjelper deplex_reads $f $raws_folder $base_path/100_formated_data/deplexed/ 7
done

samples=($(ls $base_path/100_formated_data/deplexed | cut -f1 -d"_" | sort | uniq))
