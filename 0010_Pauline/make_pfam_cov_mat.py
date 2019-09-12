import os
import pandas
from os.path import join as pjoin
from tqdm import tqdm
pfams_path = "proteom/pfams/"

def process_hmm_file(f) :
    domtblout_head = ["target_name" , "target_accession" , "query_name" , "query_accession" , "E-value","score-sequence" , "bias-sequence" , "bdE-value","score-best-domain" , "bias--best-domain" , "exp" , "reg" , "clu" , "ov" , "env" , "dom" , "rep" , "inc" , "description_of_target"]
    data = pandas.read_csv(f, delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,18))
    data = data.loc[data['bdE-value'] < 10e-9]
    pfams_dict = {p : [] for p in set(data['target_name'])}
    for a,b in data.iterrows():
        pfams_dict[b['target_name']] += [b['query_name']]
    return pfams_dict

with open("all_gffs.list") as handle :
    gff_files = handle.readlines()

prot2contig = {}
for f in tqdm(gff_files):
    with open(pjoin("../..", f[:-1])) as handle:
        for line in handle:
            if line == '##FASTA\n':
                break
            elif line[0] != "#":
                ass_type = "coass" if "1500_assemblies" in f else "single"
                assembler = "megahit" if "megahit" in f else "spades"
                cont_id = line.split("\t")[0]
                traits = {t.split("=")[0] : t.split("=")[1]  for t in line.split("\t")[-1].split(";")}
                cds_id = traits.get('ID')
                prefix = ass_type + "_" + assembler + "_"
                bad_prefix = "single_" + assembler + "_"
                cont_id  = cont_id.replace("metabat_", prefix)
                if cds_id:
                    cds_id = bad_prefix + cds_id
                    assert not prot2contig.get(cds_id)
                    prot2contig[cds_id] = cont_id

all_outs = os.listdir(pfams_path)

all_pfams = {}
for f in tqdm(all_outs):
    tt = process_hmm_file(f)
    for k, v in tt.items():
        if not all_pfams.get(k):
            all_pfams[k] = []
        all_pfams[k] += v

all_pfams_as_ctgs = {k : [ prot2contig[vv] for vv in set(v)]  for k, v in tqdm(all_pfams.items())}

covs =pandas.read_csv("cat/mapping/map_table.tsv", sep="\t")
covs.index = covs['contigName']
lengDict = covs['contigLen'].to_dict()
covs = covs[[c for c in covs.columns if c.endswith(".bam")]]

pfam_cov_dict = { k : covs.loc[v].sum() for k,v in tqdm(all_pfams_as_ctgs.items())}
pfam_pandas = pandas.DataFrame.from_dict(pfam_cov_dict, orient="index")

with open("/home/moritz/dbs/pfam/pfam_dict.csv") as handle:
    pfam2id = {l.split()[0] : l.split()[1].split(".")[0]  for l in handle}
    id2pfam = {v : k for k,v in pfam2id.items()}

with open("/home/moritz/dbs/pfam/sc_pfams.txt" ) as handle:
    sc_pfams = [id2pfam[s.strip()] for s in handle]

norm_factor = pfam_pandas.loc[sc_pfams].mean().to_dict()
normed_pfam_pandas = pandas.DataFrame.from_dict({ k : {kk : vv/norm_factor[k] for kk, vv in v.items()} for k, v in pfam_pandas.to_dict().items()})

normed_pfam_pandas.to_csv("normed_pfams.csv")

nifHs = [ "coass_megahit_FullCoassemblyNo3-1842-678" , "coass_megahit_FullCoassembly-4464-220"]
nifCovs = covs.loc[nifHs].sum()/pfam_pandas.loc[sc_pfams].mean()
