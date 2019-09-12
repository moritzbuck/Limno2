  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

load("main_var.Rdat")

station_md$station = paste("ST", station_md$station, sep="")
dna_r = fread("../map_table.tsv",h=T, sep="\t")
contig_md = dna_r[,.(contigName, contigLen)]
rm(dna_r)

rna_norm_facts = colSums(rna_reads)/mean(colSums(rna_reads))
dna_norm_facts = colSums(dna_reads)/mean(colSums(dna_reads))

bin_wise = all_table[, list(rna_logFC = median(log2FoldChange.x), dna_logFC = median(log2FoldChange.y), up_reged = sum(padj.x < 0.05 & log2FoldChange.x >0 )/length(log2FoldChange.x), down_reged = sum(padj.x < 0.05 & log2FoldChange.x <0 )/length(log2FoldChange.x)), by = bin]

bin_dat = function(bini){
  sub_info = cds_info[bin == bini & V1 %in% row.names(rna_reads), ]
  sub_counts = as.data.table(rna_reads[sub_info[,V1],])
  sub_counts$cds_id = sub_info[,V1]

  sub_data = melt(sub_counts, id = "cds_id", id_vars = "cds_id", value.name="counts", variable.name="station")
  sub_data = as.data.table(left_join(left_join(sub_data,cds_info,by= c("cds_id" = "V1")), station_md, by = c("station","station")))
  sub_data$cov = sub_data$counts / rna_norm_facts[sub_data$station] / sub_data$length

  sub_dna =   as.data.table(dna_reads[levels(as.factor(sub_data$contig)),])
  sub_dna$contig_id = levels(as.factor(sub_data$contig))

  dna_melt = melt(sub_dna, id = "contig_id", id_vars = "contig_id", value.name="dna_counts", variable.name="station")
  sub_data = as.data.table(left_join(sub_data, dna_melt, by=c("contig" = "contig_id", "station" = "station")))
  sub_data$contigLen = contig_md[match(sub_data$contig,contigName), contigLen]
  sub_data$dna_cov = sub_data$dna_counts / dna_norm_facts[sub_data$station] / sub_data$contigLen
  sub_data
}

make_cov_mat = function()
{
  dd = as.data.table(rna_reads)
  dd$cds_id = row.names(rna_reads)
  dd = melt(dd, id="cds_id", variable.name = "station",value.name="reads")
  dd = as.data.table(left_join(left_join(dd,cds_info,by= c("cds_id" = "V1")), station_md, by = c("station","station")))
  dd$cov = dd$reads / rna_norm_facts[dd$station] / dd$length
  dd_per_sample = dd[, list( median_exp = mean(cov), water.type = Water.type[1] ) , by = list(station,bin) ]
  rna_cors = sapply(levels(factor(dd_per_sample$bin)), function(x) sapply(levels(factor(dd_per_sample$bin)), function(y){
    M = dd_per_sample[bin == x]
    N = dd_per_sample[bin == y]
    M = M[match(M[,station],colnames(rna_reads))]
    N = N[match(N[,station],colnames(rna_reads))]
    cor(N$median_exp, M$median_exp)
  } ))
  dd_per_water = dd_per_sample[ , list(cov = sum(median_exp)) , by = list(bin, water.type)]
  md = data.frame(water=log10(sapply(levels(factor(dd_per_sample$bin)), function(x) dd_per_water[bin == x & water.type =="Seawater" , cov]/ dd_per_water[bin == x & water.type !="Seawater" , cov])))
  md$water[md$water == Inf] = max(md$water[md$water != Inf])
  md$water[md$water == -Inf] = min(md$water[md$water != -Inf])

  dd = as.data.table(dna_reads)
  dd$contig_id = row.names(dna_reads)
  dd = melt(dd, id="contig_id", variable.name = "station",value.name="reads")
  dd = as.data.table(left_join(left_join(dd,contig_md,by= c("contig_id" = "contigName")), station_md, by = c("station","station")))
  dd$bin = sub("-[0-9]*$","", dd$contig_id, perl = TRUE)
  dd$cov = dd$reads / dna_norm_facts[dd$station] / dd$contigLen
  dd_per_sample = dd[, list( median_exp = mean(cov), water.type = Water.type[1] ) , by = list(station,bin) ]
  dna_cors = sapply(levels(factor(dd_per_sample$bin)), function(x) sapply(levels(factor(dd_per_sample$bin)), function(y){
    M = dd_per_sample[bin == x]
    N = dd_per_sample[bin == y]
    M = M[match(M[,station],colnames(dna_reads))]
    N = N[match(N[,station],colnames(dna_reads))]
    cor(N$median_exp, M$median_exp)
  } ))
  dd_per_water = dd_per_sample[ , list(cov = sum(median_exp)) , by = list(bin, water.type)]
  md = data.frame(water=log10(sapply(levels(factor(dd_per_sample$bin)), function(x) dd_per_water[bin == x & water.type =="Seawater" , cov]/ dd_per_water[bin == x & water.type !="Seawater" , cov])))
  md$water[md$water == Inf] = max(md$water[md$water != Inf])
  md$water[md$water == -Inf] = min(md$water[md$water != -Inf])



}
