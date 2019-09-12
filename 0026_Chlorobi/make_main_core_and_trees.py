
import pandas
import os
from os.path import join as pjoin
from pynol.common.sequence.Genomic import Genomic
from pynol.common.genome.Genome import Genome
from pynol.common.taxonomy.Taxonomy import Taxonomy
from pynol.common.sequence.Feature import Feature
from tqdm import tqdm
from pynol.common.sequence.CDS import CDS
from pynol.common.sequence.RNA import RNA
from Bio import SeqIO
from pynol.common.COGs.COG import COG
from pynol.common.COGs.COGing import COGing
from mongo_thingy import connect, Thingy
from bson.objectid import ObjectId
from Bio.SeqRecord import SeqRecord
from pynol.tools.gene_prediction.Prokka import Prokka
from pynol.tools.phylogeny.Muscle import Muscle
from pynol.common.sequence.sets.Alignment import Alignment
from ete3 import Tree, NodeStyle, TreeStyle, faces
import igraph
from pynol.common.sequence.Sequence import Sequence


import uuid

connect("mongodb://localhost/chlorobi")
taxonomy = Taxonomy()
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/bac_taxonomy_r86.tsv")
genomes = list(Genome.find())
coging = COGing.find_one()

#preclustering = coging.precluster("/home/moritzbuck/people/0026_Chlorobi/all_v_all.fastani")

#main_core = coging.compute_core(cutoff = len(preclustering) / 5,precluster = "/home/moritzbuck/people/0026_Chlorobi/all_v_all.fastani")

#for c in main_core:
#    cog = COG.find_one(c)
#    cog.main_core = True
#    cog.save()

main_core = [c for c in tqdm(coging) if c.main_core]
main_core = set([c for c in main_core])
usable_core = [ c for c in  main_core if len(c._feature_list)/len(set(c.genomes)) < 1.05]
completness = lambda g : sum([g.id in c.genomes for c in useable_core])/len(useable_core)
#

#for c in tqdm(set.union(*clade_cores)):
#    for cds in tqdm(COG.find_one(c).feature_list):
#        cds.cog = c
#        cds.save()


#musc = Muscle()
#musc.make()
#musc.save()
#musc = Muscle.find_one()

#os.system("ls *.fna > genome_list.txt")


#aligns = { c : musc.run([s.sequence for s in c.feature_list], [str(s.id)for s in c.feature_list], name = str(c.id) + "_nuc") for c in usable_core}
#for f in tqdm(os.listdir("/tmp/")):
#    if f.endswith("_aligned.fasta"):
#        Alignment.FromRawFile(f.split("_")[0], "/tmp/" + f , None)


#cores_alis = [Alignment.find_one({'name' : str(l.id)})  for l in usable_core]
#cat_seqs = Alignment.cat(cores_alis, [str(g.id)for g in genomes])
#combo_ali = Alignment.FromSeqRecs(cores_alis,cat_seqs)
#combo_ali.name="concatenades_ali"
#combo_ali.save()
#combo_ali = Alignment.find_one({'name' : "concatenades_ali"})
#combo_block = combo_ali.block_seqs()
#blocked_ali = Alignment.FromSeqRecs([combo_ali],combo_block)
#blocked_ali.name = "blocked_combo_ali"
#blocked_ali.save()
blocked_ali = Alignment.find_one({'name' : "blocked_combo_ali"})

#blocked_ali.write_fasta("blocked_core.fasta")
#os.system("fasttree < blocked_core.fasta > blocked_core.fastree.tree")

core_tree =  Tree("blocked_core.fastree.tree")

get_genome = lambda s : Genome.find_one(ObjectId(Sequence.find_one( ObjectId(Sequence.find_one(ObjectId(s)).other_ids['original'])).other_ids['original']))

for l in core_tree.iter_leaves():
     genom = get_genome(l.name)
     l.features.add('species')
     l.features.add('chlorobi')
#     l.name = genom.name + "." + genom.taxonomy['gtdbtk'].replace(";",".")
     l.name = genom.id
     l.species = None


def nodes2species(node):
    genomes = [n for n in node.iter_leaf_names()]
    return set([Genome.find_one(g).species_cluster for g in genomes])

for t in core_tree.traverse():
    t.features.add('species')
    t.features.add('chlorobi')
    t.chlorobi = None
    t.species = None
    if t.name == "":
        t.species = nodes2species(t)
    else :
        t.chlorobi == Genome.find_one(ObjectId(t.name)).taxonomy['gtdbtk']

nb_chloros = lambda edge : sum(["Chlorobi" in Genome.find_one(ObjectId(t.name)).taxonomy['gtdbtk'] for t in edge.iter_leaves()])
nb_ignas = lambda edge : sum(["Ignavi" in Genome.find_one(ObjectId(t.name)).taxonomy['gtdbtk'] for t in edge.iter_leaves()])

nb_igna = nb_ignas(core_tree)
root = min([tt for tt in core_tree.traverse() if nb_ignas(tt) == nb_igna], key = lambda t : len(list(t.iter_leaves())))
core_tree.set_outgroup(root)

nb_chloro = nb_chloros(core_tree)
root = min([tt for tt in core_tree.traverse() if nb_chloros(tt) == nb_chloro], key = lambda t : len(list(t.iter_leaves())))
core_tree.set_outgroup(root)


for t in core_tree.traverse():
    t.features.add('species')
    t.features.add('chlorobi')
    t.chlorobi = None
    t.species = None
    if t.name == "":
        t.species = nodes2species(t)
    else :
        t.chlorobi == Genome.find_one(ObjectId(t.name)).taxonomy['gtdbtk']


for t in core_tree.traverse():
    if t.species:
        for tt in t.iter_descendants():
            if tt.species == t.species:
                tt.species = None



monophyl_clades_sp = { i : t.species for i, t in enumerate(core_tree.traverse()) if t.species and len(t.species) > 10}
monophyl_clades = { i : list(t.iter_leaf_names()) for i, t in enumerate(core_tree.traverse()) if t.species and len(t.species) > 10}
clade2node = { i : t for i, t in enumerate(core_tree.traverse())}
root_path = "/home/moritzbuck/people/0026_Chlorobi/"
fastani = pjoin(root_path , "all_v_all.fastani" )

clade_cores = []
for i,clade in tqdm(monophyl_clades.items()) :
    members = [Genome.find_one(g) for g in clade]
#    core = coging.compute_core(name = "clade_" + str(i) , cutoff = len(monophyl_clades_sp[i]) /10 ,precluster = fastani, genome_set=members)
    with open(pjoin(root_path, "cores", "data" , "{taxon}.core.csv".format(taxon = "clade_" + str(i)))) as handle:
        core = set([ObjectId(l[:-1]) for l in handle])
    clade_cores += [core]



for i, t in enumerate(core_tree.traverse()):
    t.features.add('full_core')
    t.features.add('ancestor_core')
    t.features.add('specific_core')
    t.features.add('id')
    t.full_core = None
    t.ancestor_core = None
    t.specific_core = None
    t.id = i

node2core = {clade2node[i] : core for i,core in zip(monophyl_clades, clade_cores)}
for n,c in node2core.items():
    n.full_core = c

main_core = node2core[core_tree]
for n in node2core:
    papis = [p.full_core for p in n.iter_ancestors() if p in node2core]
    if len(papis) > 0 :
        n.ancestor_core = set(main_core).union(set.union(*papis))
        n.specific_core = n.full_core.difference(n.ancestor_core)

    else :
        n.ancestor_core = main_core
        n.specific_core = main_core

core_tree.write(outfile="blocked_core.fastree.named.tree")


def make_tree(t, name):
    cols = {'c__Chlorobia' : "green",
             'c__Ignavibacteria' : "red",
             'c__Kapabacteria' : "orange",
             'c__UBA10030' : "purple"}

    t = t.copy()

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
#    ts.mode = "c"
    def layout(leaf):
        if not leaf.is_leaf():
            ns = NodeStyle()
            ns["hz_line_color"] = "red" if leaf.species and len(leaf.species) > 10 else "black"
            ns["hz_line_width"] = 50 if leaf.species and len(leaf.species) > 10 else 1
            ns["vt_line_width"] = len(leaf.specific_core) if leaf.specific_core else 1
            ns["size"] = 0.01 #if leaf.support > 0.95 else 0
            txt = "node_{id} ({count})".format(id = leaf.id, count = len(leaf.specific_core))  if leaf.specific_core else ""
            text = faces.TextFace(txt)
            faces.add_face_to_node(text, leaf, column=0)
            if leaf.specific_core:
                print(len(leaf.specific_core))
            leaf.set_style(ns)
        else :
            ns = NodeStyle()
            taxono = Genome.find_one(leaf.name).taxonomy['gtdbtk']
            phylum = taxono.split(";")[2]
            leaf.name = "cluster_" + str(Genome.find_one(leaf.name).species_cluster) + " " + Genome.find_one(leaf.name).name + " " + taxono
            ns['bgcolor'] = cols.get(phylum)
            ns["size"] = 0
            leaf.set_style(ns)

    ts.layout_fn = layout
    t.render(name + ".pdf", tree_style = ts)



make_tree(core_tree, pjoin(root_path, "core_tree"))

cog2function = { cog : COG.find_one(cog).get_function()[0] for cog in tqdm(set.union(*clade_cores)
)}
function2cog = {f : [c for c, ff in cog2function.items() if ff == f] for f in set(cog2function.values())}
function2cog['hypothetical protein'] = []

genome2cluster = {g.id : g.species_cluster for g in genomes}

def make_proximity_graph(temp_dict, dist_cutoff = 5000, count_cutoff = 0.5):
    all_cogs = set([f.cog for f in sum(temp_dict.values(),[])])

    all_edges = {(f,g) : 0 for f in all_cogs for g in all_cogs}
    for k, v in temp_dict.items():
        for t in v:
            for z in v:
                if t.genomic == z.genomic and abs(t.start-z.start) < dist_cutoff and t != z:
                    all_edges[(t.cog, z.cog)] += 1

    abundant_edges = {k : v for k, v in all_edges.items() if v > len(temp_dict)*count_cutoff}

    if len(abundant_edges) == 0:
        return {c : None for c in all_cogs}
    prox_graph = igraph.Graph()
    vertexDeict = { v : i for i,v in enumerate(set(sum([list(x) for x in abundant_edges.keys()], [])))}
    prox_graph.add_vertices(len(vertexDeict))
    prox_graph.add_edges([(vertexDeict[x[0]], vertexDeict[x[1]]) for x in abundant_edges.keys()])
    vertex2cog = { v : k for k,v in vertexDeict.items()}
    cog_grafs = [[vertex2cog[cc] for cc in c ] for c in prox_graph.components(mode=igraph.WEAK)]

    return {cc : str(i) if len(c) > 1 else 'None' for i,c in enumerate(cog_grafs) for cc in c }


node2cogs = {c: c.specific_core  for c in core_tree.traverse() if c.specific_core}
for taxon, core in tqdm(node2cogs.items()):
    data = {}
    genome_set = taxon.get_leaf_names()
    species_set = {genome2cluster[g]  for g in genome_set}

    cluster_set = {Genome.find_one(g).species_cluster for g in genome_set}
    cogs = [COG.find_one(ObjectId(c)) for c in core]
    temp_dict = {g : [] for g in genome_set}
    for cog in tqdm(cogs):
        funct = cog.get_function()
        feats = cog.feature_list
        count = 0
        for f in feats:
            if temp_dict.get(f.genome) != None:
                count += 1
                temp_dict[f.genome] += [f]
        found_genome = set(genome_set).intersection(cog.genomes)
        found_species = {genome2cluster[g]  for g in found_genome}
        not_found_genome = set(cog.genomes).difference(genome_set)
        not_found_species = {genome2cluster[g]  for g in not_found_genome}

        within = len(found_species)/len(species_set)
        without = len(not_found_species)/len(set(genome2cluster.values()))
#        found_fraction = len(within)/len(cluster_set)
#        found_outside = len(without)/len(genome_clusters)
        enrichment = within/without if without != 0 else float("inf")
        data[cog.id] = {
                        'name' : cog.name,
                        'function' : funct[0],
                        'fct_fract' : funct[1],
                        'hypo_fract' : funct[2],
                        'found_fraction' : within,
                        'found_outside' : without,
                        'enrichment' :  enrichment,
                        'copy_per_genome' : count/len(genome_set),
                        'functional_similar_cogs': ";".join([COG.find_one(c).name for c in function2cog[funct[0]] if c != cog.id])
                        }


    cog_clustering = make_proximity_graph(temp_dict)
    for cog, group in cog_clustering.items():
        data[cog]['operon'] = group
    pandas.DataFrame.from_dict(data, orient='index').sort_values(by="operon").to_csv(pjoin(root_path,'cores/', "node_" + str(taxon.id) + "_core_data.csv"))

genome2tax = { g.id :  g.taxonomy['gtdbtk'] for g in genomes}
for cid in tqdm(set.union(*clade_cores)):
    cog = COG.find_one(cid)
    SeqIO.write([SeqRecord(seq = c.protein, name = "", id =
"{tax}:{pid}:{fct}".format(tax = genome2tax[c.genome].replace(';','.'), pid = c.pretty_id, fct = c.more['product']), description = "") for c in cog.feature_list], pjoin(root_path, 'cogs/', "COG_" + cog.name + ".fasta") , "fasta")

genome_dat = {g.name : {'ani_cluster'  : g.species_cluster, 'taxonomy' : g.taxonomy['gtdbtk'], 'completness' : completness(g), 'assembly_length' : len(g), 'estimated_len' : len(g)/g.completness} for g in tqdm(genomes)}
pandas.DataFrame.from_dict(genome_dat, orient = 'index').sort_values(by=['taxonomy',"ani_cluster"]).to_csv(pjoin(root_path, "genome_data.csv"))

for g in tqdm(genomes):
    g.write_fasta(pjoin(root_path, "genomes", g.name + ".faa"), pretty=True)

for g in tqdm(genomes):
    g.write_fasta(pjoin(root_path, "genomes", g.name + ".fna"), pretty=True)

for g in tqdm(genomes):
    lines = [c.gff_line() for c in CDS.find({'genome' : g.id})]
    with open(pjoin(root_path, "genomes", g.name + ".gff"), "w") as handle:
        handle.writelines(lines)
