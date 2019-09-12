from Bio import SeqIO
import sys
import re

samples = {
"P6404_203" : {'fract' : 3.0, 'sample' : 905, 'water' : 'sea' },
"P6404_204" : {'fract' : 3.0, 'sample' : 909, 'water' : 'sea' },
"P6404_228" : {'fract' : 3.0, 'sample' : 910, 'water' : 'brine' },
"P6404_241" : {'fract' : 0.1, 'sample' : 898, 'water' : 'brine' },
"P6404_242" : {'fract' : 0.1, 'sample' : 902, 'water' : 'sea' },
"P6404_252" : {'fract' : 3.0, 'sample' : 911, 'water' : 'brine' },
"P6404_292" : {'fract' : 0.1, 'sample' : 912, 'water' : 'sea' },
"P6404_293" : {'fract' : 0.1, 'sample' : 916, 'water' : 'sea' }
}



def filter_fasta_by_len(args):
    length = int(args[0])
    fasta_in = args[1]
    fasta_out = args[2]

    with open(fasta_in) as handle:
        seqs = [s for s in SeqIO.parse(handle,"fasta") if len(s) > length]

    with open(fasta_out,"w") as handle:
        SeqIO.write(seqs,handle, "fasta")

def format_assembly_tables(args):
    from pandas import DataFrame

    assemstats = args[0]
    mappingstats = args[1]
    outcsv = args[2]

    with open(assemstats) as handle:
        astats =[s[:-1].split() for s in handle if s[0:8] != "filename"]

    with open(mappingstats) as handle:
        mstats = [ s[:-1].split() for s in handle]

    assemblies = set([k[0].split(".")[0] for k in astats ])
    lengths = set([k[0].split(".")[1][2:-2] for k in astats])
    lengths.add('0')
    lengths.remove('nti')

    lengths = [int(l) for l in list(lengths)]
    stats = {}
    for a in assemblies:
        for l in lengths :
            stats[(a,l)] = {}

    for a in astats:
        tt = a[0]
        ass = tt.split(".")[0]
        cutof = tt.split(".")[1][2:-2]
        cutof = int('0' if cutof == 'nti' else cutof)
        stats[(ass,cutof)]['size'] = a[1]
        stats[(ass,cutof)]['nb_contig'] = a[2]
        stats[(ass,cutof)]['min_contig_len'] = a[3]
        stats[(ass,cutof)]['median_contig_len'] = a[4]
        stats[(ass,cutof)]['mean_contig_len'] = a[5]
        stats[(ass,cutof)]['max_contig_len'] = a[6]
        stats[(ass,cutof)]['n50_len'] = a[9]
        stats[(ass,cutof)]['n90_len'] = a[11]

    for m in mstats:
        tt = m[0]
        ass = tt.split(".")[0]
        cutof = int(tt.split(".")[1][2:-2] if len(tt.split(".")) > 1 else '0')
        stats[(ass,cutof)][m[1]] = float(m[2])

    order = ['size', 'nb_contig', 'min_contig_len', 'median_contig_len', 'mean_contig_len' ,'max_contig_len','n50_len','n90_len']
    ptable = DataFrame.from_dict(stats).transpose()
    samples = [c for c in ptable.columns if c not in order]
    samples.sort()
    order = order + samples
    ptable[order].to_csv(outcsv)


def spades_assembly_qc_plot(args) :
    from pandas import DataFrame
    import pandas.rpy.common as com
    import rpy2.robjects as ro
    import rpy2.robjects.lib.ggplot2 as ggplot2

    fasta_in = args[0]
    outfile = args[1]

    get_gc = lambda s : float(s.seq.count("G")+s.seq.count("C"))/len(s.seq)
    get_cov = lambda s : float(s.description.split("_")[-1])
    get_len = lambda s : len(s.seq)

    with open(fasta_in) as handle:
        seqs = [s for s in SeqIO.parse(handle,"fasta")]

    seqs_data = DataFrame.from_dict({s.id : {'GC' : get_gc(s), 'length' : get_len(s), 'cov' : get_cov(s)} for s in seqs}).transpose()

    r_data = com.convert_to_r_dataframe(seqs_data)
    ro.r.library('ggplot2')

    x=ro.r.ggplot(r_data, ro.r.aes_string(x='length',y='cov')) + ro.r.geom_point() + ro.r.scale_y_log10() + ro.r.scale_x_log10()+ro.r.geom_vline(xintercept=1000, color="red")+ro.r.geom_vline(xintercept=2500, color="red")+ro.r.geom_vline(xintercept=10000, color="red") + ro.r.theme_bw()
    x.plot()
    ro.r('dev.copy(pdf,"%s")'%(outfile))
    ro.r('dev.off()')

def mash_dist_plot(args):
    from pandas import DataFrame
    import pandas.rpy.common as com
    import rpy2.robjects as ro
    import rpy2.robjects.lib.ggplot2 as ggplot2

    data_in = args[0]
    outfile = args[1]
    raw_data = DataFrame.from_csv(data_in)

    r_data = com.convert_to_r_dataframe(raw_data)

    ro.r.library('ggplot2')

    """
    library(vegan)

    mash_dists=as.dist(read.table("../mash_matrix.tsv"))
    mds = as.data.frame(metaMDS(mash_dists, try=100, trymax=100)$points)
    mds$assembly = row.names(mds)
    mds$assembly=sub(".contigs.fa","",mds$assembly)
    temp = as.numeric(sapply(strsplit(mds$assembly,".", fixed=TRUE), function(x) sub("bp","",sub("le","",x[2])) ))
    temp[is.na(temp)]=0
    mds$assembly = row.names(mds)
    mds$filter= as.factor(temp)
    mds$assembler = sapply(strsplit(mds$assembly,".", fixed=TRUE), function(x) if(grepl("megahit",x[1])) "megahit" else "spades" )
    mds$normed = sapply(strsplit(mds$assembly,".", fixed=TRUE), function(x) if(grepl("normed",x[1])) TRUE else FALSE )
    mds$coass = sapply(strsplit(mds$assembly,".", fixed=TRUE), function(x) if(grepl("coass",x[1])) TRUE else FALSE )
    temp=sapply(strsplit(mds$assembly,".", fixed=TRUE), function(x) tail(strsplit(x[1],"_")[[1]],1) )
    temp[temp == "coass"] = "coassembly"
    mds$sample = temp

    ggplot(mds,aes(x=MDS1, y=MDS2, col=sample, shape=assembler))+geom_point(size=4)

    pca = as.data.frame(predict(prcomp(read.tsv("distances.tsv"), sep="\t", row.names=1))[,1:4])
    pca$assembly = row.names(pca)
    pca$assembly=sub(".contigs.fa","",pca$assembly)
    temp = as.numeric(sapply(strsplit(pca$assembly,".", fixed=TRUE), function(x) sub("bp","",sub("le","",x[2])) ))
    temp[is.na(temp)]=0
    pca$assembly = row.names(pca)
    pca$filter= as.factor(temp)
    pca$assembler = sapply(strsplit(pca$assembly,".", fixed=TRUE), function(x) if(grepl("megahit",x[1])) "megahit" else "spades" )
    pca$normed = sapply(strsplit(pca$assembly,".", fixed=TRUE), function(x) if(grepl("normed",x[1])) TRUE else FALSE )
    pca$coass = sapply(strsplit(pca$assembly,".", fixed=TRUE), function(x) if(grepl("coass",x[1])) TRUE else FALSE )
    temp=sapply(strsplit(pca$assembly,".", fixed=TRUE), function(x) tail(strsplit(x[1],"_")[[1]],1) )
    temp[temp == "coass"] = "coassembly"
    pca$sample = temp
    pca = as.data.frame(predict(prcomp(read.table("../mash_matrix.tsv")))[,1:2])
    """

def rename_tree(args):
    itree = args[0]
    otree = args[1]
    if len(args) > 3 :
        taxa_file = args[3] + "/data/ppafull.tax.txt"
    else :
        taxa_file = "/home/moritz/repos/mercurial/phylophlan/data/ppafull.tax.txt"

    taxa_level = int(args[2])

    taxas = DataFrame.from_csv(taxa_file, sep="\t", header=None)
    if taxa_level > -1:
        taxas = {l[0] : l[1].split(".")[taxa_level].split("_")[-1] for l in taxas.itertuples()}
    else:
        taxas = {l[0] : "|".join([t.split("_")[-1] for t in l[1].split(".")])  for l in taxas.itertuples()}

    with open(itree) as handle:
            tree = [l for l in handle.readlines()][0]

    g_keys = [k for k in taxas.keys() if k in tree]
    for k in g_keys:
            tree=tree.replace(k,taxas[k])

    with open(otree, "w") as handle:
            handle.write(tree)



if __name__ == '__main__':
    fnct_name = sys.argv[1]
    args = sys.argv[2:]
    possibles = globals().copy()
    possibles.update(locals())
    fnct = possibles.get(fnct_name)
    if not fnct:
        raise NotImplementedError("%s does not exist" % fnct_name)
    fnct(args)
