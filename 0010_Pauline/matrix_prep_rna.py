from pandas import DataFrame, Index
from tqdm import tqdm
from Bio import SeqIO
from pandas import read_csv
from numpy import mean, median

with open("all_samples.gff") as handle:
    gff = {[vv.split("=")[1] for vv in l[:-1].split("\t")[-1].split(";") if "ID=" in vv][0] : l.split()[0]  for l in tqdm(handle) if ("all_samples-" in l.split()[0] or l.startswith("k141")) if l[0] != "#" and l[0] != ">" and "ID=" in l }

gff2 = { k : { 'contig'  : v, 'bin' : "_".join(k.split("_")[0:2]) }  for k , v in gff.items()}
gc_content = lambda s : float(str(s.seq).count('G')+str(s.seq).count('C'))/len(s)

with open("all_cdss.ffn") as handle:
    for seq in tqdm(SeqIO.parse(handle,"fasta")):
        if not gff2.get(seq.id):
            gff2[seq.id] = {}
        gff2[seq.id]['GC']=gc_content(seq)
        gff2[seq.id]['length']=len(seq)

cds_md = DataFrame.from_dict(gff2, orient="index")

fucked_up = read_csv("temp.csv")
fucked_up = fucked_up.loc[["unbinned" in v for v in fucked_up['Unnamed: 0']]]
cds_md.loc[fucked_up['Unnamed: 0'], 'contig'] = list(fucked_up['contig'])
cds_md.loc[fucked_up['Unnamed: 0'], 'bin'] = 'all_samples-unbinned'
cds_md.to_csv("cds_table.csv")
cds_md = cds_md.loc[ cds_md.contig == cds_md.contig]

lib_md = read_csv("metadata.csv")
lib_md.index = lib_md['Unnamed: 0']
lib2station = { k : str(v) for k,v in dict(lib_md['station']).items()}

dna_map = read_csv("mapping/map_table.tsv", sep="\t")
raw_dna = ((dna_map.loc[:,tuple([c for c in dna_map.columns if not "-var"  in c and "P6404" in c])]).transpose()*dna_map['contigLen']).transpose()
raw_dna = raw_dna.loc[:, tuple([ll for ll in raw_dna.columns if 'P6404_203.2' not in ll])]
raw_dna.index = dna_map['contigName']
dna_map.index = dna_map['contigName']
raw_dna.columns = Index([l.split(".")[0] for l in raw_dna.columns ])
raw_dna_per_station = raw_dna.groupby(lib2station, axis=1).sum()

normed_dna_cov = ((raw_dna_per_station/raw_dna_per_station.sum(0)).transpose()*dna_map['contigLen']).transpose()

normed_dna_fold = normed_dna_cov.apply(lambda x : x/mean(x))
normed_dna_cov.to_csv("normed_dna_cov.csv")

### Some tRNAs are missing from the cds_md, clean up eventually

rna_map = read_csv("mapping/RNA_map.tsv", sep="\t")
raw_rna = ((rna_map.loc[:,tuple([c for c in rna_map.columns if not "-var"  in c and "ST" in c])]).transpose()*rna_map['contigLen']).transpose()
raw_rna.index = rna_map['contigName'].apply(lambda s : s.split( )[0])
raw_rna.columns = Index([l.replace("_sorted.bam","").replace("ST","") for l in raw_rna.columns ])
raw_rna = raw_rna.loc[cds_md.index]

norm_facts_rns = normed_dna_fold.loc[cds_md.loc[raw_rna.index,'contig'], raw_rna.columns]
norm_facts_rns.index = raw_rna.index

norm_rna = raw_rna/norm_facts_rns
norm_rna = norm_rna.fillna(0)
norm_rna = norm_rna.loc[norm_rna.sum(1) > 0]
norm_rna.to_csv("normed_rna_counts.csv")

raw_dna_per_station.to_csv("raw_dna__per_station_counts.csv")
