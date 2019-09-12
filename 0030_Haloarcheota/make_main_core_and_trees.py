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
from Bio.Seq import Seq
import yaml

import uuid

connect("mongodb://localhost/halo")
taxonomy = Taxonomy()
#taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/gtdb_taxonomy..tsv")
taxonomy.from_gtdb_taxonomy_file("/home/moritzbuck/data/arc_taxonomy_r86.tsv")
data_folder = "/home/moritzbuck/people/0030_Haloarcheaota/"

genomes = list(Genome.find())
coging = COGing.find_one()

groups = ["p__Crenarchaeota", "p__Euryarchaeota", "p__Halobacterota"]


genomes = [g   for g in genomes if g.taxonomy['gtdb'].split(';')[1] in groups]
genome_ids = {g.id for g in genomes}

#preclustering = coging.precluster("/home/moritzbuck/people/0030_Haloarcheaota/all_v_all.archaea.fastani")

#main_core = coging.compute_core(cutoff = len(preclustering) / 5, name = "main")

#for c in main_core:
#    cog = COG.find_one(c)
#    cog.main_core = True
#    cog.save()

main_core = [c for c in tqdm(coging) if c.main_core]
main_core = set([c for c in main_core])
usable_core = [ c for c in  main_core if len(c._feature_list)/len(set(c.genomes)) < 1.05]
completness = lambda g : sum([g.id in c.genomes for c in usable_core])/len(usable_core)
#



musc = Muscle()
musc.make()
musc.save()
musc = Muscle.find_one()

#os.system("ls *.fna > genome_list.txt")


#aligns = { c : musc.run([s.sequence for s in c.feature_list], [str(s.id)for s in c.feature_list], name = str(c.id) + "_nuc") for c in usable_core}
for c in tqdm(usable_core):
     feats =  [s for s in c.feature_list if s.genome in genome_ids])
     sequences , seq_names = ([s.protein for s in feats], [str(s.id) for s in feats]
     name = str(c.id) + "_aa"
     SeqIO.write([SeqRecord(id = sname, seq = Seq(str(seq)), description="") for seq, sname in zip(sequences, seq_names)], "/tmp/" + name + ".faa", "fasta")

aligns = {c : [s for s in SeqIO.parse(pjoin(data_folder, "alis", "align_" + str(c.id) + "_aa" + ".faa"), "fasta")] for c in tqdm(usable_core)}


for f in tqdm([v for v in os.listdir(pjoin(data_folder, "alis")) if v.startswith("align_")] ):
        Alignment.FromRawFile("subaa_" + f.split("_")[1], pjoin(data_folder, "alis", f) , None)


cores_alis = [Alignment.find_one({'name' : "subaa_" + str(l.id)})  for l in usable_core]
cores_alis = [ a for a in cores_alis if a] # remove the two bugged ones
cat_seqs = Alignment.cat(cores_alis, [str(g.id)for g in genomes])
combo_ali = Alignment.FromSeqRecs(cores_alis,cat_seqs)
combo_ali.name="concatenades_ali_subaa"
combo_ali.save()
combo_ali = Alignment.find_one({'name' : "concatenades_ali_subaa"})
combo_block = combo_ali.block_seqs()
blocked_ali = Alignment.FromSeqRecs([combo_ali],combo_block)
blocked_ali.name = "blocked_combo_ali_subaa"
blocked_ali.save()
blocked_ali = Alignment.find_one({'name' : "blocked_combo_ali_subaa"})

blocked_ali.write_fasta("blocked_core.faa")
os.system("fasttree < blocked_core.faa > blocked_core.fastree.tree")

for c in tqdm(set.union(*clade_cores)):
    for cds in tqdm(COG.find_one(c).feature_list):
        cds.cog = c
        cds.save()

get_genome = lambda s : Genome.find_one(ObjectId(Sequence.find_one( ObjectId(Sequence.find_one(ObjectId(s)).other_ids['original'])).other_ids['original']))

def nodes2genomes(node):
    genomes = [n for n in node.iter_leaf_names()]
    return set(genomes)

def init_tree(tree) :
    for t in tree.traverse():
        t.features.add('genomes')
        t.features.add('taxo')
        t.taxo = None
        t.genomes = None
        if t.name == "":
            t.genomes = nodes2genomes(t)
        else :
            t.taxo = Genome.find_one(ObjectId(t.name)).taxonomy['gtdb']
    return tree



core_tree =  Tree("blocked_core.fastree.tree")
for l in core_tree.iter_leaves():
#     l.name = genom.name + "." + genom.taxonomy['gtdbtk'].replace(";",".")
     l.name = get_genome(ObjectId(l.name)).id

core_tree = init_tree(core_tree)

nb_crens = lambda edge : sum(["p__Crenarchaeota" in t.taxo for t in edge.iter_leaves()])
#nb_ignas = lambda edge : sum(["p_Halo" not in Genome.find_one(ObjectId(t.name)).taxonomy['gtdb'] for t in edge.iter_leaves()])

#nb_igna = nb_ignas(core_tree)
#root = min([tt for tt in core_tree.traverse() if nb_ignas(tt) == nb_igna], key = lambda t : len(list(t.iter_leaves())))
#core_tree.set_outgroup(root)

to_nb_cren = nb_crens(core_tree)
root = min([tt for tt in core_tree.traverse() if nb_crens(tt) == to_nb_cren], key = lambda t : len(list(t.iter_leaves())))
core_tree.set_outgroup(root)
core_tree = init_tree(core_tree)



for t in core_tree.traverse():
    if t.genomes:
        for tt in t.iter_descendants():
            if tt.genomes == t.genomes:
                tt.genomes = None



monophyl_clades_sp = { i : t.genomes for i, t in enumerate(core_tree.traverse()) if t.genomes and len(t.genomes) > 10}
monophyl_clades = { i : list(t.iter_leaf_names()) for i, t in enumerate(core_tree.traverse()) if t.genomes and len(t.genomes) > 10}
clade2node = { i : t for i, t in enumerate(core_tree.traverse())}
fastani = pjoin(data_folder , "all_v_all.archaea.fastani" )

clade_cores = []
for i,clade in tqdm(monophyl_clades.items()) :
    members = [Genome.find_one(g) for g in clade]
#    core = coging.compute_core(name = "clade_" + str(i) , cutoff = len(monophyl_clades_sp[i]) /5 ,precluster = None, genome_set=members)
    with open(pjoin(data_folder, "cores", "data" , "{taxon}.core.csv".format(taxon = "clade_" + str(i)))) as handle:
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
    cols =     {'p__' : "#5f8b74",
     'p__Altiarchaeota' : "#5e532d",
     'p__Asgardarchaeota' : "#448398",
     'p__Crenarchaeota' : "#92625d",
     'p__Euryarchaeota' : "#264133",
     'p__Hadesarchaeota' : "#d5a39d",
     'p__Halobacterota': "#3c2f1b",
     'p__JdFR-18': "#8fcfd3" ,
     'p__Micrarchaeota': "#5b2628",
     'p__Nanoarchaeota': "#c8d7b2",
     'p__Thermoplasmatota': "#311119",
     'p__UAP2': "#a4946c"}

    t = t.copy()

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
#    ts.mode = "c"
    def layout(leaf):
        if not leaf.is_leaf() or leaf.name == "":
            ns = NodeStyle()
            ns["hz_line_color"] = "red" if leaf.genomes and len(leaf.genomes) > 10 else "black"
            ns["hz_line_width"] = 50 if leaf.genomes and len(leaf.genomes) > 10 else 1
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
            taxono = leaf.taxo
            phylum = taxono.split(";")[1]
            leaf.name = "cluster_" + str(Genome.find_one(leaf.name).species_cluster) + " " + Genome.find_one(leaf.name).name + " " + taxono
            ns['bgcolor'] = cols.get(phylum, "#ffffff")
            ns["size"] = 0
            leaf.set_style(ns)

    ts.layout_fn = layout
    t.render(name + ".pdf", tree_style = ts)




make_tree(core_tree, pjoin(data_folder, "core_tree"))
clade2genome = {"clade_" + str(i) : [Genome.find_one(g).name  for g in t.genomes] for i,t in enumerate(core_tree.traverse()) if t.genomes and len(t.genomes) > 10}
with open(pjoin(data_folder, "clade2genomes.yml") , 'w') as outfile:
    yaml.dump(clade2genome, outfile, default_flow_style=False)

cog2function = { cog : COG.find_one(cog).get_function()[0] for cog in tqdm(set.union(*clade_cores)
)}
function2cog = {f : [c for c, ff in cog2function.items() if ff == f] for f in set(cog2function.values())}
function2cog['hypothetical protein'] = []

#genome2cluster = {g.id : g.species_cluster for g in genomes}

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

with open("/home/moritzbuck/data/ec2kegg/kegg_maps.txt") as handle:
    kegg2name = {l.split()[0] : " ".join(l.split()[1:]) for l in  handle.readlines() if l.startswith("    ")}

with open("/home/moritzbuck/data/ec2kegg/kegg_path2ec.txt") as handle:
    ec2path = [(" ".join(l.split()[1:]).replace("ec:","") , l.split()[0].replace("path:ec","")) for l in  handle.readlines()]

ec2path = {ec : [tt[1] for tt in ec2path if tt[0] == ec] for ec in set([t[0] for t in ec2path])}
ec2path = {ec : [kegg2name[vv] for vv in v] for ec, v in ec2path.items() }

node2cogs = {c: c.specific_core  for c in core_tree.traverse() if c.specific_core}
for taxon, core in tqdm(node2cogs.items()):
    data = {}
    genome_set = taxon.get_leaf_names()

    cogs = [COG.find_one(ObjectId(c)) for c in core]
    temp_dict = {g : [] for g in genome_set}

    for cog in tqdm(cogs):
        funct = cog.get_function()
        ecs = cog.ec_number()

        feats = cog.feature_list
        count = 0
        for f in feats:
            if temp_dict.get(f.genome) != None:
                count += 1
                temp_dict[f.genome] += [f]
        found_genome = set(genome_set).intersection(cog.genomes)
        not_found_genome = set(cog.genomes).difference(genome_set)

        within = len(found_genome)/len(genome_set)
        without = len(not_found_genome)/(len(genomes)-len(genome_set)) if (len(genomes)-len(genome_set)) > 0 else float("inf")
#        found_fraction = len(within)/len(cluster_set)
#        found_outside = len(without)/len(genome_clusters)
        enrichment = within/without if without != 0 else float("inf")
        data[cog.id] = {
                        'name' : cog.name,
                        'function' : funct[0],
                        'fct_fract' : funct[1],
                        'hypo_fract' : funct[2],
                        'ec_number' : ecs[0].replace("hypothetical protein",""),
                        'ec_fract' : ecs[1],
                        'noec_fract' : ecs[2],
                        'pathways' : ";".join(ec2path.get(ecs[0], [""])),
                        'found_fraction' : within,
                        'found_outside' : without,
                        'enrichment' :  enrichment,
                        'copy_per_genome' : count/len(genome_set),
                        'functional_similar_cogs': ";".join([COG.find_one(c).name for c in function2cog[funct[0]] if c != cog.id])
                        }


    cog_clustering = make_proximity_graph(temp_dict)
    for cog, group in cog_clustering.items():
        data[cog]['operon'] = group
    pandas.DataFrame.from_dict(data, orient='index').sort_values(by="operon").to_csv(pjoin(data_folder,'cores/', "node_" + str(taxon.id) + "_core_data.csv"))

genome2tax = { g.id :  g.taxonomy['gtdb'] for g in Genome.find()}
for cid in tqdm(set.union(*clade_cores)):
    cog = COG.find_one(cid)
    SeqIO.write([SeqRecord(seq = c.protein, name = "", id =
"{tax}:{pid}:{fct}".format(tax = genome2tax[c.genome].replace(';','.'), pid = c.pretty_id, fct = c.more['product']), description = "") for c in cog.feature_list], pjoin(data_folder, 'cogs/', "COG_" + cog.name + ".fasta") , "fasta")

def clade_faa(node ):
    genomes = set(node.genomes)
    cogs = node.specific_core
    seqs = []
    for cid in cogs:
        cog = COG.find_one(cid)
        seqs += [SeqRecord(seq = c.protein, name = "", id = "{tax}:{pid}:{fct}".format(tax = genome2tax[c.genome].replace(';','.'), pid = c.pretty_id, fct = c.more['product']), description = "") for c in cog.feature_list if c.genome in genomes]
    return seqs

for i, t in enumerate(core_tree.traverse()):
    if t.specific_core:
        SeqIO.write(clade_faa(t), pjoin(data_folder, 'cores/node_{i}.faa'.format(i=i) ) ,"fasta")

genome_dat = {g.name : {'ani_cluster'  : g.species_cluster, 'taxonomy' : g.taxonomy['gtdb'], 'completness' : completness(g), 'assembly_length' : len(g), 'estimated_len' : len(g)/completness(g)} for g in tqdm(genomes)}
pandas.DataFrame.from_dict(genome_dat, orient = 'index').sort_values(by=['taxonomy',"ani_cluster"]).to_csv(pjoin(data_folder, "genome_data.csv"))

for g in tqdm(genomes):
    g.write_fasta(pjoin(data_folder, "genomes", g.name + ".faa"), pretty=True)

for g in tqdm(genomes):
    g.write_fasta(pjoin(data_folder, "genomes", g.name + ".fna"), pretty=True)

for g in tqdm(genomes):
    lines = [c.gff_line() for c in CDS.find({'genome' : g.id})]
    with open(pjoin(data_folder, "genomes", g.name + ".gff"), "w") as handle:
        handle.writelines(lines)
