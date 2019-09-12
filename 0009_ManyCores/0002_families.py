import os, sys
from os.path import join as pjoin
from tqdm  import tqdm
from joblib import Parallel, delayed
from pandas import DataFrame
import random

home = os.environ['HOME']
sys.path.append(pjoin(home,"repos/moritz/MiComPy/"))

import micompy.databases.refseq.refseq
from micompy.common.tools.workbench import WorkBench


def single_ANI(genome_tuple):
    g1 = genome_tuple[0]
    g2 = genome_tuple[1]
    return ((g1.name,g2.name) , g1.get_ANI(g2))


bench = WorkBench()
bench.default_bench()

refseq_db = micompy.databases.refseq.refseq.Refseq(workbench = bench)
num_cores = 12

all_complete = [g for g in tqdm(refseq_db) if  g.metadata['assembly_level'] == 'Complete_Genome' and len(g.pfams) > 0]

all_families = set([g.get_taxo().get('family') for g in all_complete])

fam_clusters = {f : [g for g in all_complete if g.get_taxo().get('family') == f] for f in tqdm(all_families)}
big_families = {f : v for f,v in fam_clusters.items() if len(v) > 50}
subset_big_fams = {f : random.sample(v,50) for f,v in big_families.items()}
all_gs = sum(subset_big_fams.values(),[])

pfam_sim = lambda p,q : float(len(p.pfams.intersection(q.pfams)))/max(len(p.pfams),len(q.pfams)) if max(len(p.pfams),len(q.pfams)) != 0 else None
pfam_simi = {}


for c in tqdm(subset_big_fams.values()):
    for g1 in tqdm(c):
        for g2 in c:
            if not pfam_simi.has_key((g1,g2)):
                pfam_simi[(g1,g2)] = pfam_sim(g1,g2)

for g in tqdm(all_gs):
    if not os.path.exists(pjoin(g.path, "ref")):
        bench.tools['BBMap'].make_index(g)

DataFrame.from_dict({(k[0].name, k[1].name) : {'pfam_simi' : v} for k,v in pfam_simi.items()}).transpose().to_csv("pfam_simis.csv")
ANIs = {}
tt = DataFrame.from_csv("ANIs_fams.csv")
ANIs = {(t[0],t[1][0]) : {'ANI' : t[1][1], 'coverage' : t[1][2]} for t in tt.iterrows()}


for c in tqdm(subset_big_fams.values()):
    to_compute = set()
    for g1 in tqdm(c):
        for g2 in c:
            if not (g1,g2) in to_compute and not (g2,g1) in to_compute and not (g1.name,g2.name) in ANIs.keys() and not (g2.name,g1.name) in ANIs.keys():
                to_compute.add((g1,g2))
    data = Parallel(n_jobs=num_cores)(delayed(single_ANI)(i) for i in tqdm(to_compute))
    ANIs.update(data)
    DataFrame.from_dict(ANIs).transpose().to_csv("ANIs_fams.csv")
