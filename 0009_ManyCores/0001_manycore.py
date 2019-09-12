import os, sys
from os.path import join as pjoin
from tqdm  import tqdm
from joblib import Parallel, delayed
from pandas import DataFrame

home = os.environ['HOME']
sys.path.append(pjoin(home,"repos/moritz/MiComPy/"))

import micompy.databases.refseq.refseq
from micompy.common.tools.workbench import WorkBench

def compu_pfam(g):
    g.get_pfams()

def pre_process():
    to_process = [g for g in refseq_db if g.metadata['assembly_level'] == "Complete_Genome"]
    Parallel(n_jobs=num_cores)(delayed(compu_pfam)(i) for i in tqdm(to_process))

def parse_mcl_clusters(file):
    with open(file) as handle:
        clusters = {l[:-1].split()[0] : l[:-1].split() for l in handle.readlines()}
        clusters = {k : [refseq_db[vv] for vv in v] for k,v in clusters.items()}
        return clusters

def single_ANI(genome_tuple):
    g1 = genome_tuple[0]
    g2 = genome_tuple[1]
    return ((g1.name,g2.name) , g1.get_ANI(g2))

bench = WorkBench()
bench.default_bench()

refseq_db = micompy.databases.refseq.refseq.Refseq(workbench = bench)
num_cores = 12

chosen_clusters = parse_mcl_clusters(pjoin(refseq_db.metadata_path,"mcl_genome_clusters","mcl_clusters_" + str(80) + "_I_" + str(1.4) + ".tsv"))
chosen_clusters = {k : [vv for vv in v if vv] for k,v in chosen_clusters.items()}

pfam_sim = lambda p,q : float(len(p.pfams.intersection(q.pfams)))/max(len(p.pfams),len(q.pfams)) if max(len(p.pfams),len(q.pfams)) != 0 else None
pfam_simi = {}

for c in tqdm(chosen_clusters.values()):
    for g1 in tqdm(c):
        for g2 in c:
            if not pfam_simi.has_key((g1,g2)) and not pfam_simi.has_key((g2,g1)):
                pfam_simi[(g1,g2)] = pfam_sim(g1,g2)

ANIs = {}
to_compute = set()
for c in tqdm(chosen_clusters.values()):
    for g1 in tqdm(c):
        for g2 in c:
            if not (g1,g2) in to_compute and not (g2,g1) in to_compute:
                to_compute.add((g1,g2))
