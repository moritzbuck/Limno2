library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
options(width=200))

md = read.csv("../metadata.csv", h=TRUE, row.names=1, as.is=TRUE)
row.names(md) = sub("P6404_","",row.names(md))

station_md = read.csv("../RNA_station_dat.csv", sep=",", as.is=TRUE)
to_clean = names(which(apply(station_md, 2, function(x) sum(grepl("±", x))) > 0))
for(c in to_clean){
  station_md[,c] = sapply(strsplit(station_md[,c]," "), function(x) as.numeric(x[1]))
}
station_md$station = sapply(strsplit(station_md$Station.nr," "),tail,1)
station_md = as.data.table(station_md)

rna_raw = fread("../RNA_map.tsv",h=T, sep="\t")
rna_raw$annot = sapply(strsplit(rna_raw$contigName," "),function(x) paste(x[2:length(x)], collapse=" ", sep=" ") )
rna_raw$cds_id = sapply(strsplit(rna_raw$contigName," "),head, 1 )
cds_md = rna_raw[,.(cds_id, contigLen,annot)]
rna_raw = rna_raw[,c("cds_id",grep("var", grep("ST", colnames(rna_raw),v=TRUE), invert=TRUE, value=TRUE)), with=FALSE]



melted_rna = melt(rna_raw, value.name="coverage", variable.name = "sample", id.vars="cds_id")
melted_rna$station = sub("_sorted.bam","",sub("ST","",melted_rna$sample))
melted_rna = merge(melted_rna, cds_md, by = "cds_id")

melted_rna$tot_bases = melted_rna$coverage*melted_rna$contigLen
sample_rna_base_count = melted_rna[,.(total_bases = sum(tot_bases)), by =station]
norm_fact = sample_rna_base_count[,.(norm_fact=total_bases/mean(total_bases))][,norm_fact]
names(norm_fact) = sample_rna_base_count[,station]

melted_rna$normed_bases = melted_rna$tot_bases / norm_fact[melted_rna$station]
melted_rna$normed_coverage = melted_rna$normed_bases / melted_rna$contigLen

melted_rna$bin = sapply(strsplit(melted_rna$cds_id,"_"), function(x) paste(head(x,2), sep="_", collapse="_"))
melted_rna$bin[grep("k",melted_rna$bin)] = "unbinned"

dna_raw = fread("../map_table.tsv",h=T, sep="\t")
contig_md = dna_raw[,.(contigName, contigLen)]
dna_raw = dna_raw[,c("contigName",grep("var", grep("P64", colnames(dna_raw),v=TRUE), invert=TRUE, value=TRUE)), with=FALSE]
null_samples = names(which (colSums(dna_raw[,grep("P64",colnames(dna_raw), val=TRUE), with = FALSE]) == 0))
dna_raw = dna_raw[, !colnames(dna_raw) %in% null_samples, with = FALSE]

melted_dna = melt(dna_raw, value.name="coverage", variable.name = "sample", id.vars="contigName")
melted_dna$library = sub(".bam","",sub("P6404_","",melted_dna$sample))
melted_dna = merge(melted_dna, contig_md, by = "contigName")

melted_dna$tot_bases = melted_dna$coverage*melted_dna$contigLen
sample_dna_base_count = melted_dna[,.(total_bases = sum(tot_bases)), by =library]
dna_norm_fact = sample_dna_base_count[,.(norm_fact=total_bases/mean(total_bases))][,norm_fact]
names(dna_norm_fact) = sample_dna_base_count[,library]

melted_dna$normed_bases = melted_dna$tot_bases / dna_norm_fact[melted_dna$library]
melted_dna$normed_coverage = melted_dna$normed_bases / melted_dna$contigLen
melted_dna = melted_dna[library != "203.2"]

melted_dna$bin = sapply(strsplit(melted_dna$contigName,"-"), function(x) paste(head(x,2), sep="-", collapse="-"))
melted_dna$bin[grep("k",melted_dna$bin)] = "unbinned"
melted_dna$station = md[melted_dna$library,"station"]
melted_dna$type = md[melted_dna$library,"type"]

dna_per_sample = summarize(group_by(melted_dna,contigName, station, type), total_coverage = sum(normed_coverage))
dna_per_sample = as.data.table(dna_per_sample)
dna_per_sample$bin = sapply(strsplit(dna_per_sample$contigName,"-"), function(x) paste(head(x,2), sep="-", collapse="-"))
dna_per_sample$bin[grep("k",dna_per_sample$bin)] = "unbinned"
dna_per_sample$station  = as.character(dna_per_sample$station)
dna_per_sample = merge(dna_per_sample,station_md[,.(LightNitrogen,DarkNitrogen, station, LightC, DarkC)], by="station")

dna_per_bin = summarize(group_by(melted_dna,library, bin), coverage = sum(normed_bases)/sum(contigLen))
dna_per_bin_per_site = summarize(group_by(merge(dna_per_bin, md, by="library", by.y="row.names"),station,type, bin), coverage=sum(coverage))
dna_per_bin_per_site$station = as.character(dna_per_bin_per_site$station)

dna_nitr_cor = summarize(group_by(dna_per_sample, contigName),
  light_N_cor = cor(total_coverage,LightNitrogen, method = "spearman"),
  dark_N_cor = cor(total_coverage, DarkNitrogen, method = "spearman"),
  light_C_cor = cor(total_coverage,LightC, method = "spearman"),
  dark_C_cor = cor(total_coverage,DarkC, method = "spearman"),
  bin = bin[1]
)
dna_nitr_cor = as.data.table(dna_nitr_cor)
dna_nitr_cor = dna_nitr_cor[!is.na(light_C_cor)]

tt = summarize(group_by(dna_nitr_cor, bin), med = mean(light_N_cor))
dna_med_cors = tt$med
names(dna_med_cors) = tt$bin
dna_nitr_cor_filt = dna_nitr_cor[dna_nitr_cor$bin %in% names(dna_med_cors)]
dna_nitr_cor_filt$bin = factor(dna_nitr_cor_filt$bin, levels=names(sort(dna_med_cors)))


rna_per_sample = merge(melted_rna,station_md[,.(LightNitrogen,DarkNitrogen, station, LightC, DarkC)], by="station")
rna_per_sample_with = left_join(rna_per_sample, dna_per_bin_per_site, by = c("station" = "station", "bin" = "bin") )
rna_per_sample_with = filter(rna_per_sample_with, !is.na(coverage.y))
rna_per_sample_with = transform(rna_per_sample_with, rel.coverage = normed_coverage/coverage.y)

rna_nitr_cor = summarize(group_by(rna_per_sample_with, cds_id),
  light_N_cor = cor(normed_coverage,LightNitrogen, method = "spearman", use = "complete.obs"),
  dark_N_cor = cor(normed_coverage, DarkNitrogen, method = "spearman", use = "complete.obs"),
  light_C_cor = cor(normed_coverage,LightC, method = "spearman", use = "complete.obs"),
  dark_C_cor = cor(normed_coverage,DarkC, method = "spearman", use = "complete.obs"),
  light_N_rel = cor(rel.coverage,LightNitrogen, method = "spearman", use = "complete.obs"),
  dark_N_re = cor(rel.coverage, DarkNitrogen, method = "spearman", use = "complete.obs"),
  light_C_rel = cor(rel.coverage,LightC, method = "spearman", use = "complete.obs"),
  dark_C_rel = cor(rel.coverage,DarkC, method = "spearman", use = "complete.obs"),
  light_N_dna = cor(coverage.y,LightNitrogen, method = "spearman", use = "complete.obs"),
  dark_N_dna = cor(coverage.y, DarkNitrogen, method = "spearman", use = "complete.obs"),
  light_C_dna = cor(coverage.y,LightC, method = "spearman", use = "complete.obs"),
  dark_C_dna = cor(coverage.y,DarkC, method = "spearman", use = "complete.obs"),
  mean_rel_cov = mean(rel.coverage),
  mean_dna_cov = mean(coverage.y),
  mean_rna_cov = mean(normed_coverage),
  bin = bin[1]
)

rna_nitr_cor = as.data.table(rna_nitr_cor)
rna_nitr_cor = rna_nitr_cor[!is.na(light_N_cor)]

tt = summarize(group_by(rna_nitr_cor, bin), med = mean(light_N_rel))
rna_med_cors = tt$med
names(rna_med_cors) = tt$bin
rna_nitr_cor_filt = rna_nitr_cor[rna_nitr_cor$bin %in% names(rna_med_cors)]
rna_nitr_cor_filt$bin = factor(rna_nitr_cor_filt$bin, levels=names(sort(rna_med_cors)))
tata = merge(as.data.frame(dna_med_cors), as.data.frame(rna_med_cors), by="row.names")

ggplot( tata, aes(x=dna_med_cors, y=rna_med_cors, label=Row.names))+geom_point() +geom_label_repel(data=filter(tata,Row.names == shewi | rna_med_cors >0.4))



cds_info = fread("../cds_table.csv")
cds_md = merge(cds_md, cds_info, by = "cds_id", by.y = "X")
bin_md = fread("../bin_stats.csv", sep=",", head=T)
bin_md$avg_prot_len = bin_md$length*bin_md$coding_density/bin_md$nb_proteins
