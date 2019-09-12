import os
import gzip
from tqdm import tqdm
from Bio.Seq import Seq
import pandas

os.chdir("/home/moritz/people/0006_Erken_Genome_Project/10000_SAGs/")
undet_folder =  "../0000_rawdata/0100_SAG_data/0110_libtest/Undertermined/"
os.listdir(undet_folder)

lib_one = [f for f in os.listdir(undet_folder) if f.endswith(".fastq.gz") if "L001" in f and "I1" in f][0]
lib_two = [f for f in os.listdir(undet_folder) if f.endswith(".fastq.gz") if "L002" in f and "I1" in f][0]

dif = lambda s1, s2 : sum([ a != b for a,b in zip(s1,s2)])
hh = lambda v : v[1]

def most_sim(seq, seq_set, max_err = 3):
    difs = [(qes,dif(seq, qes)) for qes in seq_set]
    vals = min(difs, key = hh)
    return vals if vals[1] < max_err else None

with gzip.open(undet_folder + "/" + lib_one) as handle:
    idxs_one = [l.decode().split(":")[-1][:-1] for l in tqdm(handle) if l.decode().startswith("@ST-")]

idxs_one = [ tuple(i.replace("N", "").split("+")) for i in idxs_one]

with open("metad.csv") as handle:
    lines = [l for l in handle.readlines() if "NXT" in l]

fwd_idxs = { l.split()[0] : l.split()[-2] for l in lines[1:]}
rev_idxs = { l.split()[0] : l.split()[-1] for l in lines[1:]}
rev_fwd_idxs = {k : [v for v, kk in fwd_idxs.items() if kk == k] for k in set(fwd_idxs.values())}
rev_rev_idxs = {k : [v for v, kk in rev_idxs.items() if kk == k] for k in set(rev_idxs.values())}
fwd_idxs_list = list(fwd_idxs.values())
rev_idxs_list = list(rev_idxs.values())


uniq_idxs = set(idxs_one)

fwd_corrected = {idx : most_sim(idx[0], fwd_idxs_list) for idx in tqdm(uniq_idxs)}
rev_corrected = {idx : most_sim(idx[1], rev_idxs_list) for idx in tqdm(uniq_idxs)}


counts = { i : 0 for i in tqdm(uniq_idxs)}
for i in tqdm(idxs_one):
    counts[i]  = counts[i] + 1

undets = {s.split("_")[0] : 0 for s in all_samples}
undets['undet'] = 0
for u in uniq_idxs:
    sami = 1
    fwd_c = fwd_corrected[u]
    rev_c = rev_corrected[u]
    if fwd_c and rev_c:
        sami = 2
        undets['undet'] += counts[u]
    if fwd_c :
        fwd_samples = rev_fwd_idxs[fwd_c[0]]
        for s in fwd_samples:
            undets[s] += counts[u]/len(fwd_samples)/sami
    if rev_c :
        rev_samples = rev_rev_idxs[rev_c[0]]
        for s in rev_samples:
            undets[s] += counts[u]/len(rev_samples)/sami

undet_dict = {}
undet_dict['counts'] = undets
undet_dict['cell'] = {k : k.split('-')[2] if k != 'undet' else k for k in undet_dict['counts'].keys()}
undet_dict['kit'] = {k : k.split('-')[-1] if k != 'undet' else k for k in undet_dict['counts'].keys()}
pandas.DataFrame.from_dict(undet_dict).to_csv("undet_data_NXT.csv")






with open("metad.csv") as handle:
    lines = [l for l in handle.readlines() if "NXT" not in l]

fwd_idxs = { l.split()[0] : l.split()[-2] for l in lines[1:]}
rev_idxs = { l.split()[0] : l.split()[-1] for l in lines[1:]}
rev_fwd_idxs = {k : [v for v, kk in fwd_idxs.items() if kk == k] for k in set(fwd_idxs.values())}
rev_rev_idxs = {k : [v for v, kk in rev_idxs.items() if kk == k] for k in set(rev_idxs.values())}
fwd_idxs_list = list(fwd_idxs.values())
rev_idxs_list = list(rev_idxs.values())


with gzip.open(undet_folder + "/" + lib_two) as handle:
    idxs_two = [l.decode().split(":")[-1][:-1] for l in tqdm(handle) if l.decode().startswith("@ST-")]

idxs_two = [ tuple(i.replace("N", "").split("+")) for i in idxs_two]
uniq_idxs = set(idxs_two)

fwd_corrected = {idx : most_sim(idx[0], fwd_idxs_list) for idx in tqdm(uniq_idxs)}
rev_corrected = {idx : most_sim(idx[1], rev_idxs_list) for idx in tqdm(uniq_idxs)}


counts = { i : 0 for i in tqdm(uniq_idxs)}
for i in tqdm(idxs_two):
    counts[i]  = counts[i] + 1

undets = {s.split("_")[0] : 0 for s in all_samples}
undets['undet'] = 0
for u in uniq_idxs:
    sami = 1
    fwd_c = fwd_corrected[u]
    rev_c = rev_corrected[u]
    if fwd_c and rev_c:
        sami = 2
        undets['undet'] += counts[u]
    if fwd_c :
        fwd_samples = rev_fwd_idxs[fwd_c[0]]
        for s in fwd_samples:
            undets[s] += counts[u]/len(fwd_samples)/sami
    if rev_c :
        rev_samples = rev_rev_idxs[rev_c[0]]
        for s in rev_samples:
            undets[s] += counts[u]/len(rev_samples)/sami

undet_dict = {}
undet_dict['counts'] = undets
undet_dict['cell'] = {k : k.split('-')[2] if k != 'undet' else k for k in undet_dict['counts'].keys()}
undet_dict['kit'] = {k : k.split('-')[-1] if k != 'undet' else k for k in undet_dict['counts'].keys()}
pandas.DataFrame.from_dict(undet_dict).to_csv("undet_data.csv")












all_reads = os.listdir("1000_reads_proc/1200_trimmomatic/")
all_samples = [r.split(".")[0].replace("_1P", "") for r in all_reads if r.endswith("_1P.fastq.gz")]

all_asses = [a.replace(".fasta", "") for a in os.listdir("4000_assembly_analysis/4999_genomes/") if a.endswith("_100cells_2500bp.fasta")]

def get_stats(fasta, sample):
    dups_template = "4000_assembly_analysis/4300_mapping/{fasta}/bams/{sample}_sorted.stats".format(sample = sample, fasta = fasta)
    nodups_template = "4000_assembly_analysis/4300_mapping/{fasta}/bams/{sample}.stats".format(sample = sample, fasta = fasta)
    with open(dups_template) as handle:
        dups_vals = [float(l[:-1].split("(")[-1].split("%")[0]) for l in handle if "%" in l]
    with open(nodups_template) as handle:
        mapped, proper, singles = tuple([float(l[:-1].split("(")[-1].split("%")[0]) for l in handle if "%" in l])
    with open(dups_template) as handle:
        reads = [float(l[:-1].split()[0]) for l in handle if "paired in"in l][0]
    with open(nodups_template) as handle:
        dupt_fract = 1-[float(l[:-1].split()[0])/reads for l in handle if "paired in"in l][0]

    out = { 'map' : mapped/100,
            'proper_map' : proper/100,
            'unpaired' : singles/100,
            "dups" : dupt_fract,
            'cell' : fasta.split("-")[2],
            'kit' : fasta.split("_")[0].split("-")[-1] ,
            'read_cell' : sample.split("-")[2],
            'read_kit' : sample.split("_")[0].split("-")[-1],
            }

    out['same_cell'] = out['cell'] == out['read_cell']
    out['same_kit'] = out['kit'] == out['read_kit']


    return out

map_stats = { (f, s) : get_stats(sample = s, fasta = f) for f in all_asses for s in all_samples}
pandas.DataFrame.from_dict(map_stats, orient="index").to_csv("map_stats.csv")
