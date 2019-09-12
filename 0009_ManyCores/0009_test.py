import os, sys
from os.path import join as pjoin
from tqdm  import tqdm
home = os.environ['HOME']
sys.path.append(pjoin(home,"repos/moritz/MiComPy/"))
from micompy.databases.database import Database
from micompy.common.genome import Genome
from micompy.common.tools.workbench import WorkBench

bench = WorkBench()
bench.default_bench()
g = Genome("GCF_000005845.2", ".", "GCF_000005845.2.fna", workbench = bench)
Db = Database("test_db", [g])
Db.process()
test = bench['HMMer'].hmmsearch_pfam_presence(g)
