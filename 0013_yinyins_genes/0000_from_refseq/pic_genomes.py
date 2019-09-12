import os, sys
from os.path import join as pjoin
from tqdm  import tqdm
from joblib import Parallel, delayed
from pandas import DataFrame

home = os.environ['HOME']
sys.path.append(pjoin(home,"repos/moritz/MiComPy/"))

import micompy.databases.refseq.refseq
from micompy.common.tools.workbench import WorkBench

bench = WorkBench()
bench.default_bench()

refseq_db = micompy.databases.refseq.refseq.Refseq(workbench = bench)


with open("/home/moritz/taxas.csv") as handle:
    taxas = [" ".join(t.split()[:2]) for t in handle.readlines()]

taxas = set(taxas)

sub_refseq = [g for g in refseq_db if " ".join(g.metadata.get('organism_name').split()[:2]) in taxas]


picked = set([" ".join(g.metadata.get('organism_name').split()[:2]) for g in sub_refseq])
taxas.difference(picked)

base = "/media/moritz/LENOVO_USB_HDD/genomes_for_rnf/"

for g in tqdm(sub_refseq):
    if not os.path.exists(pjoin(base," ".join(g.metadata.get('organism_name').split()[:2]))):
        os.makedirs(pjoin(base," ".join(g.metadata.get('organism_name').split()[:2])))
    shutil.copy(pjoin(g.path, "original_files" , g.metadata.get('ftp_path').split("/")[-1] + ".fna.gz" ), pjoin(base," ".join(g.metadata.get('organism_name').split()[:2]))) 
