import sys
import os
from os.path import join as pjoin
from pandas import DataFrame
from pandas import Index
from tqdm  import tqdm
from numpy import nan
from pylab import *
from numpy import sort

home = os.environ['HOME']
sys.path.append(pjoin(home,"repos/moritz/MiComPy/"))

from micompy.pipes.analyses import *
from micompy.common.genome import Genome
from micompy.common.utils.renaming_tree import renaming_tree
from micompy.common.utils.iotl_annotations import *


root = pjoin(home, "Data/people/0008_bacties")
data_root = pjoin(root, "000_data/")
pre_anal_root = pjoin(root, "1000_preanalysis")
checkm_dir = pjoin(pre_anal_root, "1100_checkm/")
if not os.path.exists(checkm_dir):
    os.makedirs(checkm_dir)


cpus = 11

all_closed = os.listdir(pjoin(data_root, "closed"))
all_cultfree = os.listdir(pjoin(data_root, "cultivation_frees"))

metadata = DataFrame.from_csv(pjoin(data_root,"taxontable18602_09-jun-2017.xls"), sep = "\t")
metadata.index = Index(metadata.index.to_series().apply(str))
metadata['short_name'] = nan #metadata['taxon_oid'].apply(str)
# check if all genomes are in the metadata sheet
assert all([g[:-6] in metadata.index for g in all_closed + all_cultfree])

metadata = metadata.transpose().to_dict()

all_closed_genomes = [ Genome(g[:-6], pjoin(data_root,"processed_genomes",g[:-6]), pjoin(data_root,'closed', g ), metadata[g[:-6]]) for g in all_closed]
all_cultfree_genomes = [ Genome(g[:-6], pjoin(data_root,"processed_genomes",g[:-6]), pjoin(data_root,'cultivation_frees', g ), metadata[g[:-6]]) for g in all_cultfree]

all_genomes = all_closed_genomes + all_cultfree_genomes

def run():
    annotation(all_genomes, cpus=cpus)
    checkm(all_genomes, cpus=cpus, output = pjoin(checkm_dir, "checkm_all.txt"))
    phylophlan(all_genomes, cpus=cpus, output = pjoin(checkm_dir, "phyluo.txt"))
    renaming_tree(pjoin(checkm_dir, "phyluo.txt"), pjoin(checkm_dir, "phyluo_renamed.txt"), {  g.name : "|".join([g.metadata['Class'],g.metadata['Order'],g.metadata['Family'], g.metadata['Genus']] )for g in all_genomes })renaming_tree(pjoin(checkm_dir, "phyluo.txt"), pjoin(checkm_dir, "phyluo_renamed.txt"), {  g.name : "|".join([g.metadata['Class'],g.metadata['Order'],g.metadata['Family'], g.metadata['Genus']] )for g in all_genomes })

    clade_data(all_genomes, {g.name : g.metadata['Class'] for g in all_genomes}, pjoin(checkm_dir, "phylophlan_class.txt")  )
    class_data(all_genomes, {g.name : g.metadata['Culture Type'] for g in all_genomes if g.metadata['Culture Type'] == g.metadata['Culture Type']}, pjoin(checkm_dir, "phylophlan_type.txt")  )
