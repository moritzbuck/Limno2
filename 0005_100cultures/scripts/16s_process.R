library(plyr)
library(ggplot2)
library(vegan)
library(pheatmap)
library(igraph)
library(RColorBrewer)


parse_taxon_str = function(x)
{
  x=gsub(" ","_",x)
  x=gsub(".","_",x,fixed=T)
  x=gsub("-","_",x,fixed=T)
  x=gsub("[","_",x,fixed=T)
  x=gsub("]","_",x,fixed=T)

  val = sub("[A-Za-z0-9:_/()]*[(]([0-9_]*)[)]","\\1", x, perl=TRUE)
  val= as.numeric(gsub("_",".",val,fixed=T))

  string = sub("([A-Za-z:_/0-9]*)[(][0-9_]*[)]","\\1", x, perl=TRUE)
  names(val) = string
  val
}

get_taxa =  function(otus, taxa, level = NA, cutof=0.6, na.val="tip")
{
  tt = taxa[otus, "V2"]

  taxa_list = strsplit(gsub("\\\"", "" ,tt), ",")
  parsed_taxa_list = lapply(taxa_list, parse_taxon_str)
  parsed_taxa_list = lapply(parsed_taxa_list, function(x) x[!grepl("uncultured",names(x))])
  parsed_taxa_list = lapply(parsed_taxa_list, function(x) x[!grepl("Unknown",names(x))])
  final_taxa_list = lapply(parsed_taxa_list, function(x) names(x > cutof))

  if(is.na(level))
    output = sapply(final_taxa_list,tail,1)
  else
  {
    if(na.val=="tip")
      output = sapply(final_taxa_list,function(x)  if(length(x) < (level+1)) tail(x,1) else x[level] )
    else
      output = sapply(final_taxa_list,function(x)  if(length(x) < (level+1)) na.val else x[level] )
  }
  names(output) = otus
  output
}


collapse_otu_table = function(table, taxa, level)
{
  tax = get_taxa(row.names(table), taxa, level)
  odat = ddply(cbind(as.data.frame(table),tax),.(tax), colwise(sum))
  row.names(odat) = odat[,1]
  odat = odat[,-1]
  odat
}
per_otu_cutoff = 50
per_sample_cutoff = 300

data = read.table("jamt_otu_table_97_nosingletons.txt",h=T, comment.char="&", row.names=1, sep="\t")

full.data.jamt = data
jamt.taxa = read.table("jamt_97_singletons.rdp.sintax", row.names=1, as.is=T,sep="\t")
#taxa = read.table("jamt_unoise.sintax", row.names=1, as.is=T,sep="\t")
row.names(jamt.taxa) = sapply(strsplit(row.names(jamt.taxa),";"),head,1)

negs = grep("NTC",colnames(data),v=T)
if(length(negs) == 1){
  suspicious_otus = row.names(data)[data[,negs] > per_otu_cutoff]
} else {
  suspicious_otus = row.names(data)[rowSums(data[,negs]) > per_otu_cutoff]
}

data[data < per_otu_cutoff] = 0
data = data[, colSums(data) > per_sample_cutoff]
data = data[rowSums(data) > 0, ]
freq_dat = t((t(data)/colSums(data)))

jamt = freq_dat[,grepl("Jamtl",colnames(freq_dat))]
jamt = jamt[rowSums(jamt)>0,]
jamt = jamt[!row.names(jamt) %in% suspicious_otus,]

ana = freq_dat[,grepl("Ana",colnames(freq_dat))]
ana = ana[rowSums(ana)>0,]
ana = ana[!row.names(ana) %in% suspicious_otus,]


data = read.table("erken_otu_table_97_nosingletons.txt",h=T, comment.char="&", row.names=1, sep="\t")

full.data.erken = data

erken.taxa = read.table("erken_97_singletons.rdp.sintax", row.names=1, as.is=T,sep="\t")
row.names(erken.taxa) = sapply(strsplit(row.names(erken.taxa),";"),head,1)

negs = grep("NTC",colnames(data),v=T)
if(length(negs) == 1){
  suspicious_otus = row.names(data)[data[,negs] > per_otu_cutoff]
} else {
  suspicious_otus = row.names(data)[rowSums(data[,negs]) > per_otu_cutoff]
}

data[data < per_otu_cutoff] = 0
data = data[, colSums(data) > per_sample_cutoff]
data = data[rowSums(data) > 0, ]
freq_dat = t((t(data)/colSums(data)))

erken = freq_dat[,grepl("erken",colnames(freq_dat))]
erken = erken[rowSums(erken)>0,]
erken = erken[!row.names(erken) %in% suspicious_otus,]




phmap = function(table,taxa, level = NA, bad_otus = NULL){
  if(is.na(level))
  {
    f_table = table[!row.names(table) %in% bad_otus,]
    row.names(f_table) = paste(row.names(f_table),get_taxa(row.names(f_table), taxa, NA))
  }
  else{
    level = if(level=="tips") NA else level
    f_table = collapse_otu_table(table = table[!row.names(table) %in% bad_otus,], taxa,level)
  }

  solos = as.vector(apply(f_table[, which(colSums(f_table>0) == 1 )],2,function(x) names(which(x > 0))))

  adj = apply(f_table > 0, 1, function(x) apply(f_table > 0, 1, function(y) sum(x & y)/sum(x)))
  net = graph_from_adjacency_matrix(adj, weighted = "freq", diag = FALSE)

  return(list(table = f_table, net=net, solos = solos ))
}

make_heatmap = function(table, metadat = NULL)
{
  pheatmap((table > 0)+0 , clustering_distance_rows="binary",
  clustering_distance_cols="binary", treeheight_row=4, treeheight_col=4, legend=FALSE,annotation_legend = FALSE,
  color = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100), border_color=NA, annotation_row=metadat, show_colnames =FALSE)
}

make_network = function(obj, cutoff = 0.8,bottom.cutoff = 1, top.cutoff = 0.6, interactive = FALSE, layout=layout_nicely, pval=FALSE, pval.cutoff=0.001)
{
  net = obj$net
  table = obj$table

  tt =as_edgelist(net)
  pvals = sapply(1:nrow(tt), function(x) edge.prob(table, tt[x,1], tt[x,2], 10))
#    pvals = p.adjust(pvals, method="fdr")
  E(net)$pvals = pvals

  if(pval)
  {

    net = delete.edges(net, which(E(net)$pvals > pval.cutoff))
  }
  else
  {
  net = delete.vertices(net,which(rowSums(table[names(V(net)),]> 0) < bottom.cutoff )  )
  net = delete.edges(net, which(E(net)$freq < cutoff))
  net = delete.vertices(net,which(rowMeans(table[names(V(net)),]> 0) > top.cutoff))
  }

  comps = components(net)
  connected = names(comps$membership[which(comps$membership %in% which(comps$csize != 1))])

  net = induced.subgraph(net,connected)

  sizes = 10*log10(rowSums(table[names(V(net)),]> 0))
  names(sizes) = names(V(net))
  width=(E(net)$freq - cutoff)*9 +1
  colors = rep(NA, length(V(net)))
  colors[which(names(V(net)) %in% obj$solos)] = "orange"
  edge.colors = rep("gray", length(E(net)))

  edge.colors[E(net)$pvals < pval.cutoff] = "red"

  if(interactive)
  {
      tkplot(net, vertex.size= sizes,vertex.label.color = "black", edge.width = width, layout=layout)
  }
  else
  {
    plot(net, vertex.color= colors, vertex.size= sizes,vertex.label.color = "black", edge.color= edge.colors, edge.width = width, layout = layout)
  }
  return(net)
}


load(file="erken.coords")
erken.adj = phmap(table=erken, erken.taxa, level=5)
erken.net = make_network(erken.adj, bottom.cutoff=3, top.cutoff=0.5, cutoff=0.5, layout=erken.coords)
dev.copy(svg,"Erken_network.svg", width = 12, height=12)
dev.off()
erken.meta = data.frame(row.names = row.names(erken.adj$table), kept = row.names(erken.adj$table) %in% names(V(erken.net))+0, isolates = row.names(erken.adj$table) %in% erken.adj$solos+0 )
make_heatmap(erken.adj$table, metadat = erken.meta)
dev.copy(svg,"Erken_heatmap.svg", width = 12, height=12)
dev.off()
dev.off()


load(file="jamt.coords")
jamt.adj = phmap(table=jamt, jamt.taxa, level=5)
jamt.net = make_network(jamt.adj, bottom.cutoff=3, top.cutoff=0.5, cutoff=0.5, layout=jamt.coords)
dev.copy(svg,"Jamt_network.svg", width = 12, height=12)
dev.off()
jamt.meta = data.frame(row.names = row.names(jamt.adj$table), kept = row.names(jamt.adj$table) %in% names(V(jamt.net))+0, isolates = row.names(jamt.adj$table) %in% jamt.adj$solos+0 )
make_heatmap(jamt.adj$table, metadat = jamt.meta)
dev.copy(svg,"Jamt_heatmap.svg", width = 12, height=12)
dev.off()
dev.off()


load(file="ana.coords")
ana.adj = phmap(table=ana, jamt.taxa, level=5)
ana.net = make_network(ana.adj, bottom.cutoff=3, top.cutoff=0.5, cutoff=0.5, layout=ana.coords)
dev.copy(svg,"Anaeroic_network.svg", width = 12, height=12)
dev.off()
ana.meta = data.frame(row.names = row.names(ana.adj$table), kept = row.names(ana.adj$table) %in% names(V(ana.net))+0, isolates = row.names(ana.adj$table) %in% ana.adj$solos+0)
make_heatmap(ana.adj$table, metadat = ana.meta)
dev.copy(svg,"Anaerobic_heatmap.svg", width = 12, height=12)
dev.off()


edge.prob = function(freqs, node1, node2, ncells){
  prob.dist = rowMeans(freqs)
  sub.bool = (freqs > 0)[c(node1,node2),]
  prob1 = (1-(1-prob.dist[node1])^ncells) # prob of picking at least one node1
  prob2 = (1-(1-prob.dist[node2])^ncells)
  prob.co = prob1*prob2 #### a bit wrong
  prob.1only = prob1*(1-prob2)
  prob.2only = prob2*(1-prob1)
  total.rolls = ncol(freqs)
  count.co = sum(colSums(sub.bool) == 2)
  count.1only = sum(sub.bool[1,] > sub.bool[2,])
  count.2only = sum(sub.bool[2,] > sub.bool[1,])
  prob.of.event = 1-(1-prob.co)^(total.rolls - count.co)
  return(prob.of.event)
}





pretty_freq_dat = freq_dat[!row.names(freq_dat) %in% suspicious_otus,]
row.names(pretty_freq_dat) = paste(row.names(pretty_freq_dat),get_taxa(row.names(pretty_freq_dat), 6))




adjacency = apply(pretty_freq_dat > 0.0, 1, function(x) apply(pretty_freq_dat > 0, 1, function(y) sum(x & y)/mean(sum(x), sum(y))) )

f_table = collapse_otu_table(table = freq_dat[!row.names(freq_dat) %in% suspicious_otus,],4)
adjacency = apply(f_table > 0.0, 1, function(x) apply(f_table > 0, 1, function(y) sum(x & y)/mean(sum(x), sum(y))) )

pheatmap(((f_table > 0) +0), fontsize=6, clustering_distance_rows="binary",
  clustering_distance_cols="binary", treeheight_row=4, treeheight_col=4, legend=FALSE,
  color = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100), border_color=FALSE)#,
  cellwidth=7, cellheight=7, width=10, height=3, filename="Erken_cultures.pdf")

pheatmap(((pretty_freq_dat > 0) +0), fontsize=7, clustering_distance_rows="binary",
  clustering_distance_cols="binary", treeheight_row=4, treeheight_col=4, legend=FALSE,
  color = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100), border_color=FALSE)#,
  cellwidth=7, cellheight=7, width=12, height=6, filename="Erken_cultures_otus.pdf")
