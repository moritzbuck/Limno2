import os


checkm_fields = ['Bin Id',
    'Marker lineage',
    '# genomes',
    '# markers',
    '# marker sets',
    '0',
    '1',
    '2',
    '3',
    '4',
    '5+',
    'Completeness',
    'Contamination',
    'Strain heterogeneity']


checkm_file = "good_bins_from_singles.checkm"
with open(checkm_file) as handle:
    dat = handle.readlines()

bin_md = {l.split()[0] : {k : v for k,v in zip(checkm_fields, l.replace(" (UID", "_(UID").split())} for l in dat[3:-1]}

sample_metadata = "/home/moritz/people/0010_Pauline/metadata.csv"
with open(sample_metadata) as handle:
    dat = handle.readlines()

sample_md = {l.split(",")[0] : {k: v for k,v in zip(["fraction", "sample", "biome", "site"],l[:-1].split(",")[1:])} for l in dat}

bins = [f[:-4] for f in os.listdir(".") if f[-4:] == ".fna"]

fulll_md = {b :  {** sample_md[b.split("-")[0]], **bin_md[b], **{'bin_id' : b.split("-")[1]} } for b in bins}

trans_table = [ "{bin}\t{biome}-{site}-{fract}-{bin_id}---Complete-{compl}-Cont-{cont}\n".format(bin=k, fract = v['fraction'], biome = v['biome'], site=v['site'], bin_id = v['bin_id'], compl=v['Completeness'], cont = v['Contamination']) for k,v in fulll_md.items()]

out_file = "good_bins_from_singles.md"
with open(out_file,"w") as handle:
     handle.writelines(trans_table)
