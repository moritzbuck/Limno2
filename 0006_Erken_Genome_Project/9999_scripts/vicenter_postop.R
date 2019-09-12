library(phyloseq)
library(ggplot2)
library(vegan)
library(scales)
library(grid)
library(metagenomeSeq)
library(randomForest)
library(magrittr)
library(Biostrings)
require(tidyr)
require(reshape2)
require(Rmisc)
library(RDPutils  )
require(seqinr)
library(DESeq2)
library(plyr)
library(dplyr)
library(qualpalr)
library(data.table)


tax <- read.csv(file="dada2_silva.taxonomy", header = TRUE)
rownames(tax) <- tax[,1]
tax <- tax[,-1]
tax <- as.matrix(tax)

otu <- read.table("otu_table.txt", sep="\t", header=TRUE, row.names=1)
otu <- otu_table(otu, taxa_are_rows = TRUE, errorIfNULL = TRUE)

mapfile = "sample_map_ErkVic.csv"
map <- read.csv(mapfile, sep=";")
rownames(map) <- map[,1]
map$RO = map$organics=="Reverse Osmosis"  | map$organics == "Mixed" | map$organics == "DOC-Shading"
map$H = map$organics=="Humin Feed"  | map$organics == "Mixed"
map$Shading =  map$organics == "Shading" | map$organics == "DOC-Shading"
map$dark = NA
map$dark[map$Shading | map$H] = "darkened"
map$dark[!(map$Shading | map$H)] = "light"
map$dark = factor(map$dark, levels = c('light', 'darkened') )
map$DOC = NA
map$DOC[map$RO] = "extra DOC"
map$DOC[!map$RO] = "no addition"
map$DOC = factor(map$DOC, levels = c('no addition', 'extra DOC') )

map <- sample_data(map[,-1])

physeq <- phyloseq(otu,tax_table = tax_table(tax))
physeq <- merge_phyloseq(physeq,map)
#physeq <- subset_taxa(physeq, Order!="Chloroplast")
physeq = subset_samples(physeq, sample_type != "negative" & sample_type != "mock")
physeq = subset_samples(physeq, experiment != "Spring" & experiment != "20 mesocosms")
physeq = subset_samples(physeq, original_name !="ExpB 1.11")

tt = tax_table(physeq)[,"Order"] != "Chloroplast"
physeq = prune_taxa( row.names(tt)[which(tt)], physeq)



# obj = subset_samples(physeq, experiment =="Experiment A"); fract= 0.01; level= 5
tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus" )

level = 4
fract = 0.01
rel_abs = physeq %>%  tax_glom(taxrank = tax_levels[level]) %>% transform_sample_counts(function(x) {x/sum(x)} )
keepers = names(which(apply(otu_table(rel_abs), 1, max) > fract))
kept = prune_taxa(keepers,rel_abs)
melted = psmelt(kept)
# transforms matrix into one-line-per-observation matrix, ggplot loves that

melted$full_tax = apply(melted[,tax_levels[1:level]], 1 ,paste, sep=";", collapse=";")
temp_map = as.vector(melted$Class[!duplicated(melted$full_tax)])
names(temp_map) = melted$full_tax[!duplicated(melted$full_tax)]

melted$full_tax = factor(melted$full_tax, levels=c(names(temp_map)[order(temp_map)], "others"))
groupings = sapply(levels(as.factor(melted$Class)), function(x) levels(factor(melted$full_tax[melted$Class == x])))

start = qualpal(length(levels(melted$Class)), list(h=c(0,360), s = c(0.7,0.7), l=c(0.6, 0.6)))$HSL[,1]
names(start) = levels(melted$Class)
cols = lapply(names(start), function(x) if(length(groupings[[x]]) > 1) qualpal(n = length(groupings[[x]]), list(h = c(start[x], start[x]), s = c(0.6, 0.9), l = c(0.6, 0.85)))$hex else hsv(start[[x]]/360, 0.7,0.6) )
cols = c(unlist(cols))
greys = "#888888"
cols = c( cols, greys)
names(cols) = c(as.vector(unlist(groupings)), "others")


make_bar_stack = function(obj, fract, level)
{
  rel_abs = obj %>%  tax_glom(taxrank = tax_levels[level]) %>% transform_sample_counts(function(x) {x/sum(x)} )
  keepers = names(which(apply(otu_table(rel_abs), 1, max) > fract))
  kept = prune_taxa(keepers,rel_abs)
  melted = psmelt(kept)
  # transforms matrix into one-line-per-observation matrix, ggplot loves that

  melted$full_tax = apply(melted[,tax_levels[1:level]], 1 ,paste, sep=";", collapse=";")
  temp_map = as.vector(melted$Class[!duplicated(melted$full_tax)])
  names(temp_map) = melted$full_tax[!duplicated(melted$full_tax)]

  melted$full_tax = factor(melted$full_tax, levels=c(names(temp_map)[order(temp_map)], "others"))

  melted = merge(melted,ddply(melted, .(DOC,dark,replicates),  summarize, Abundance = 1-sum(Abundance), full_tax = "others"), all=TRUE)

  cols = cols[levels(factor(melted$full_tax))]

  legend_text = levels(melted$full_tax)
  legend_text = sapply(strsplit(legend_text, ";"), tail, 1)

  pp = ggplot(melted, aes(x = factor(replicates), y = Abundance, fill = full_tax, order = Class)) + geom_bar(stat = "identity")+scale_fill_manual(values = cols, labels = legend_text) + facet_grid(DOC~dark)+ylim(0,1)+theme_bw()+guides(fill=guide_legend(ncol=1,name = "taxonomical group"))
  list(plot = pp, data = melted)
}

bmft = function(obj){

  deseq_res = phyloseq_to_deseq2(obj, ~dark*DOC)
  deseq_res = DESeq(deseq_res, test="Wald", fitType="parametric")
  contrasts = resultsNames(deseq_res)[-1]


  all = lapply(contrasts, function(cont){
    res =  results(deseq_res, cooksCutoff = FALSE, contrast=list(cont))
    res$tax = apply(tax_table(obj), 1 ,paste, sep=";", collapse = ";")
    res$tax = paste(res$tax,row.names(res), sep=";")
    res$level = "OTU"
    res = as.data.table( res)
    bla =
    lapply(2:6, function(x){
    tdat = obj %>%  tax_glom(taxrank = tax_levels[x])
    deseq_res = phyloseq_to_deseq2(tdat, ~dark*DOC)
    deseq_res = DESeq(deseq_res, test="Wald", fitType="parametric")
    res =  results(deseq_res, cooksCutoff = FALSE, contrast=list(cont))
    res$tax = apply(tax_table(tdat), 1 ,paste, sep=";", collapse = ";")
    res$level = tax_levels[x]
    as.data.table(res)
    }
    )

    bla = c(list(res) , bla)
    big_Tab = do.call("rbind", bla)
    big_Tab$contrast=cont
    big_Tab
  }
  )
  big_Tab = do.call("rbind", all)
  big_Tab$padj=p.adjust(big_Tab$pvalue)
  big_Tab
}


taxo_plot = function(obj, level, taxon)
{
  tt = obj %>% tax_glom(taxrank = level) %>% transform_sample_counts(function(x) {x/sum(x)} )
  keep = tax_table(tt)[,level] == taxon
  keep = row.names(keep)[which(keep)]
  tt = prune_taxa(keep,tt)
  melted = psmelt(tt)
  ggplot(melted, aes(x=organics, y=Abundance, col=organics))+geom_boxplot()+geom_jitter(width=0.5)+scale_y_log10()+ggtitle(taxon)
}

hm_data = function(obj, level = NA, tax_filt = NULL)
{
  print(level)
  mat = obj %>% filter_taxa(function(x) sum(x) > 0 , prune = TRUE)
  if(!is.na(level) & level != "OTU")
  {
      mat = mat %>% tax_glom(taxrank = level)
  } else
  {
    level ="OTU"
  }
  tt = mat %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt
  if(!is.null(tax_filt))
  {
    tt = tt[tt[,level] %in% tax_filt,]
  }
  mat = mat %>% rarefy_even_depth
  mat2 = t(apply(otu_table(mat), 1, function(x) x/mean(x)))
  if(!is.na(level) & level != "OTU")  row.names(mat2) = as.vector(tax_table(mat)[row.names(mat2),level])
  mat2 = mat2[row.names(mat2) %in% tax_filt,]

  row.order <- hclust(dist(log(mat2+1)))$order # clustering
  row.order = rev(rownames(mat2)[row.order])

  col.order <- hclust(dist(log(t(mat2+1))))$order
  col.order = colnames(mat2)[col.order]
  tt$level = level
  tt$taxon = as.vector(tt[,level])
  tt2 = ddply(tt, .(taxon), mutate, abundance = Abundance/mean(Abundance))

  tt2$Sample = factor(tt2$Sample, level = col.order)
  tt2$taxon = factor(tt2$taxon, level = row.order)

  tt2 = tt2[, !colnames(tt2) %in% tax_levels]
  tt2
}


dsplit = function(tax_str){
  vecti = strsplit(tax_str,";")[[1]]
  vecti[vecti != "NA"]
}

expA = subset_samples(physeq, experiment =="Experiment A" & !is.na(replicates) & time_point =="t5")
extt = make_bar_stack(obj = expA, fract= fract, level= level)
extt[['plot']]+ggtitle("Experiment A")+xlab("replicates")
ggsave("barstack_expA.pdf", width = 18, height = 12)

bla = bmft(expA %>% rarefy_even_depth)
sigs.expA = bla[padj < 0.001][order(padj)]

good_tax = lapply(c(tax_levels, "OTU"), function(x) sapply(strsplit(sigs.expA[level == x]$tax , ";"), function(y) tail(y[y!="NA"], 1)))
names(good_tax) = c(tax_levels, "OTU")
expA.hm = lapply(c(tax_levels[2:6], "OTU"), function(x) hm_data(expA, x, tax_filt = good_tax[[x]]))
tt = do.call(rbind,expA.hm)
levels(tt$taxon) = unlist(lapply(expA.hm, function(x) levels(x$taxon)))
tt$level = factor(tt$level, levels = c(tax_levels, "OTU"))
tt$organics = factor(tt$organics,levels(tt$organics)[c(1,4,2,3)])
ggplot(tt, aes(x=Sample, y = taxon, fill=abundance)) +geom_tile()+facet_grid(level~organics, scale="free")+scale_fill_gradient2(trans="log2", high = "blue4", low = "khaki", midpoint = 0.25,  na.value ="lightgrey")+  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave("heatmap_expA.pdf", width = 8, height = 20)

fwrite(sigs.expA, "significant_expA.csv")
sigs.expA.nonOTUs = sigs.expA[level != "OTU"]
single_plots = lapply(1:nrow(sigs.expA.nonOTUs), function(n){gg=taxo_plot(expA, sigs.expA.nonOTUs[n]$level, tail(dsplit(sigs.expA.nonOTUs[n]$tax),1) ); print(n); gg })

pdf("all_plots_expA.pdf")
for(p in single_plots)
{
  plot(p)
}
dev.off()

expB = subset_samples(physeq, experiment =="Experiment B" & !is.na(replicates) & time_point =="t5")
extt = make_bar_stack(obj = expB, fract= fract, level= level)
extt[['plot']]+ggtitle("Experiment B")+xlab("replicates")
tdat = ddply(extt[['data']], .(Order, interaction(DOC,dark)), summarize, mean_ab = mean(Abundance), sd_ab=sd(Abundance))
colnames(tdat)[2] = "fact"
ggsave("barstack_expB.pdf", width = 18, height = 12)
bla = bmft(expB %>% rarefy_even_depth)
sigs.expB = bla[padj < 0.001][order(padj)]

good_tax = lapply(c(tax_levels, "OTU"), function(x) sapply(strsplit(sigs.expB[level == x]$tax , ";"), function(y) tail(y[y!="NA"], 1)))
names(good_tax) = c(tax_levels, "OTU")
expB.hm = lapply(c(tax_levels[2:6], "OTU"), function(x) hm_data(expB, x, tax_filt = good_tax[[x]]))
tt = do.call(rbind,expB.hm)
levels(tt$taxon) = unlist(lapply(expB.hm, function(x) levels(x$taxon)))
tt$level = factor(tt$level, levels = c(tax_levels, "OTU"))
tt$organics = factor(tt$organics,levels(tt$organics)[c(1,3,4,2)])
ggplot(tt, aes(x=Sample, y = taxon, fill=abundance)) +geom_tile()+facet_grid(level~organics, scale="free")+scale_fill_gradient2(trans="log2", high = "blue4", low = "khaki", midpoint = 0.250,  na.value ="lightgrey")+  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave("heatmap_expB.pdf", width = 8, height = 20)

fwrite(sigs.expB, "significant_expB.csv")

sigs.expB.nonOTUs = sigs.expB[level != "OTU"]
single_plots = lapply(1:nrow(sigs.expB.nonOTUs), function(n){gg=taxo_plot(expB, sigs.expB.nonOTUs[n]$level, tail(dsplit(sigs.expB.nonOTUs[n]$tax),1) ); print(n); gg })

#pdf("all_plots_expB.pdf")
#for(p in single_plots)
#{
#  plot(p)
#}
#dev.off()
