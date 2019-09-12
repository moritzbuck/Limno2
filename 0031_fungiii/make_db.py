from ete3 import NCBItaxa

ncbi = NCBITaxa()

with open("metadat.csv") as handle:
    name2file = { l.split("\t")[0] : l.split("\t")[1][:-1]  for l in handle}

species = list(set([n.split()[0] for n in name2file.keys()]))


spec2tax = {k : dict([ (ncbi.get_rank([t])[t], ncbi.get_taxid_translator([t])[t]) for t in ncbi.get_lineage(v[0])]) for k,v in ncbi.get_name_translator(species).items()}

cols = c("#4bc3b7","#c949a7","#56c468","#9061cd","#85b937","#6272c0","#b9ae39","#6ca2da","#ca512f","#408f66","#d84165","#478833","#d484be","#a2b069","#a34c6c","#d6903b","#776e29","#c87a60", "grey", "darkgrey")
