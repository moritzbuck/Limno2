  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library('BiocParallel')
  library(doMC)

register(MulticoreParam(30),default=TRUE)

md = read.csv("../metadata.csv", h=TRUE, row.names=1, as.is=TRUE)

station_md = read.csv("../RNA_station_dat.csv", sep=",", as.is=TRUE)
to_clean = names(which(apply(station_md, 2, function(x) sum(grepl("±", x))) > 0))
for(c in to_clean){
  station_md[,c] = sapply(strsplit(station_md[,c]," "), function(x) as.numeric(x[1]))
}
station_md$station = sapply(strsplit(station_md$Station.nr," "),tail,1)
row.names(station_md) =  sub("Station ", "ST", station_md$Station.nr)

rna_reads = read.csv("normed_rna_counts.csv", row.names=1, as.is=TRUE)
colnames(rna_reads) = sub("X","ST",colnames(rna_reads))

rna_reads = rna_reads[,row.names(station_md)]
rna_reads = replace(round(rna_reads), rna_reads == Inf, 0)
rna_reads = rna_reads[rowSums(rna_reads) != 0,]

weirdo = "all_samples-unbinned_345951"
rna_reads = rna_reads[row.names(rna_reads) != weirdo,]

rna_raw = fread("../RNA_map.tsv",h=T, sep="\t")
rna_raw$annot = sapply(strsplit(rna_raw$contigName," "),function(x) paste(x[2:length(x)], collapse=" ", sep=" "))
cds_md = rna_raw[,.(contigName, contigLen,annot)]
cds_md$cds_id = sapply(strsplit(cds_md$contigName," "),function(x) x[1])
cds_md = data.frame(cds_md,row.names=cds_md$cds_id)
cds_md = cds_md[row.names(rna_reads),]
cds_md$bin = sapply(strsplit(cds_md$cds_id,"_"), function(x) paste(head(x,2), sep="_", collapse="_"))
cds_md$bin[grep("k",cds_md$bin)] = "unbinned"
rm(rna_raw)

bin_md = fread("../bin_stats.csv", sep=",", head=T)


cds_md = as.data.table(cds_md)
cds_info = fread("../cds_table.csv",h=T, sep=",")
cds_info = left_join(cds_info,cds_md[match(cds_id, cds_info$V1), .(cds_id,annot)], by = c("V1" = "cds_id"))
cds_info = as.data.table(cds_info)

station_md$bad_lat = as.numeric(sub("'", "", sub("º ",".", sub(".","", sub("N ", "", station_md$Latitude), fixed =T))))
#seqed = DESeqDataSetFromMatrix(countData=rna_reads, colData=station_md, design=~Salinity+bad_lat+LightNitrogen+DarkNitrogen)
seqed = DESeqDataSetFromMatrix(countData=rna_reads, colData=station_md, design=~LightNitrogen)



dd_rna = DESeq(seqed,  parallel=TRUE)
dd_rna_shrink <- lfcShrink(dd_rna, coef=2, type="ash", parallel=TRUE)
rna_table_light = as.data.table(dd_rna_shrink)
rna_table_light$ids = row.names(dd_rna_shrink)
rna_table_light = left_join(rna_table_light, cds_info, by = c("ids" = "V1"))

bla = apply(rna_reads, 1, function(x) cor(x,station_md$LightNitrogen))
rna_cor = as.data.table(bla)
colnames(rna_cor) = "log2FoldChange"
rna_cor$ids = names(bla)
rna_cor$padj=1

bla = apply(dna_reads, 1, function(x) cor(x,station_md$LightNitrogen))
dna_cor = as.data.table(bla)
colnames(dna_cor) = "log2FoldChange"
dna_cor$ids = names(bla)
dna_cor$padj=1

rna_cor = left_join(rna_cor, cds_info, by = c("ids" = "V1"))
all_table_cor = left_join(rna_cor, dna_cor, by = c("contig" = "ids"))
all_table_cor = as.data.table(all_table_cor)


### dna analysis

dna_reads = read.csv("raw_dna__per_station_counts.csv", row.names=1, as.is=TRUE)
colnames(dna_reads) = sub("X","ST",colnames(dna_reads))

dna_reads = dna_reads[,row.names(station_md)]
dna_reads = replace(round(dna_reads), dna_reads == Inf, 0)
dna_reads = dna_reads[rowSums(dna_reads) != 0,]

#seqed_dna = DESeqDataSetFromMatrix(countData=dna_reads, colData=station_md, design=~Salinity+bad_lat+LightNitrogen+DarkNitrogen)
seqed_dna = DESeqDataSetFromMatrix(countData=dna_reads, colData=station_md, design=~LightNitrogen)


dd_dna = DESeq(seqed_dna,  parallel=TRUE)

temp_dna =  lfcShrink(dd_dna, coef=2, type="ash", parallel=TRUE)
dna_table_light = as.data.table(temp_dna)
dna_table_light$ids = row.names(temp_dna)
all_table = left_join(rna_table_light, dna_table_light, by = c("contig" = "ids"))
all_table = as.data.table(all_table)

sub_main = all_table[sample(.N,floor(nrow(all_table)/10))]
main_plot =  ggplot(sub_main, aes(y=log2FoldChange.x, x=log2FoldChange.y))+geom_density2d(bins=200)+xlim(c(-max(sub_main[,log2FoldChange.y])),max(sub_main[,log2FoldChange.y]))+ylim(c(-max(sub_main[,log2FoldChange.x])),max(sub_main[,log2FoldChange.x]))

sub_main_cor = all_table_cor[sample(.N,floor(nrow(all_table_cor)/10))]

main_plot_cor =  ggplot(sub_main_cor, aes(y=log2FoldChange.x, x=log2FoldChange.y))+geom_density2d(bins=50)

pmain_plot =  ggplot(sub_main, aes(y=-sign(log2FoldChange.x)*log10(pvalue.x), x=-sign(log2FoldChange.y)*log10(pvalue.y)))+geom_density2d(bins=50)
bin_plot = function(bini)
{
  info = paste(c(colnames(bin_md), "expressed ratio"), c(bin_md[V1 == bini], nrow(all_table[bin == bini]) / nrow(cds_info[bin == bini])), sep = ": ", collapse = "\n")


  sub = filter(all_table, bin == bini)
  sub$sig = (sub$padj.x < 0.05 | sub$padj.y < 0.05)
#ggplot(sub_main, aes(y=log2FoldChange.x, x=log2FoldChange.y))+geom_point(size=0.1, col="black")
  main_plot + geom_point(data = sub, aes(col=as.factor(sig)), size =4)+
  geom_vline(xintercept=median(sub$log2FoldChange.y),col="red")+geom_hline(yintercept=median(sub$log2FoldChange.x),col="red")+geom_vline(xintercept=0, col="black")+geom_hline(yintercept=0, col="black")+
  theme_bw()+ggtitle(bini)+scale_color_manual(values=c("orange1","indianred"),na.value="gray")+xlab("genomic fold change (log2)")+ylab("transcriptomic fold change (log2)")+geom_text(data=as.data.frame(info), x=-1.5,y=2, aes(label = info))

}

pbin_plot = function(bini)
{
  info = paste(c(colnames(bin_md), "expressed ratio"), c(bin_md[V1 == bini], nrow(all_table[bin == bini]) / nrow(cds_info[bin == bini])), sep = ": ", collapse = "\n")

  sub = filter(all_table, bin == bini)
  sub$sig = (sub$padj.x < 0.05 | sub$padj.y < 0.05)
#ggplot(sub_main, aes(y=log2FoldChange.x, x=log2FoldChange.y))+geom_point(size=0.1, col="black")
  pmain_plot + geom_point(data = sub, aes(col=as.factor(sig)), size =4)+
  geom_vline(xintercept=median(-sign(sub$log2FoldChange.y)*log10(sub$pvalue.y)),col="red")+geom_hline(yintercept=median(-sign(sub$log2FoldChange.x)*log10(sub$pvalue.x)),col="red")+geom_vline(xintercept=0, col="black")+geom_hline(yintercept=0, col="black")+
  theme_bw()+ggtitle(bini)+scale_color_manual(values=c("orange1","indianred"),na.value="gray")+xlab("genomic log10(pval)")+ylab("transcriptomic log10(pval)")+geom_text(data=as.data.frame(info), x=1.5,y=3, aes(label = info))

}

cor_plot = function(bini)
{
  info = paste(c(colnames(bin_md), "expressed ratio"), c(bin_md[V1 == bini], nrow(all_table_cor[bin == bini]) / nrow(cds_info[bin == bini])), sep = ": ", collapse = "\n")


  sub = filter(all_table_cor, bin == bini)
  sub$sig = (sub$padj.x < 0.05 | sub$padj.y < 0.05)
#ggplot(sub_main, aes(y=log2FoldChange.x, x=log2FoldChange.y))+geom_point(size=0.1, col="black")
  main_plot_cor + geom_point(data = sub, aes(col=as.factor(sig)), size =4)+
  geom_vline(xintercept=median(sub$log2FoldChange.y),col="red")+geom_hline(yintercept=median(sub$log2FoldChange.x),col="red")+geom_vline(xintercept=0, col="black")+geom_hline(yintercept=0, col="black")+
  theme_bw()+ggtitle(bini)+scale_color_manual(values=c("orange1","indianred"),na.value="gray")+xlab("genomic correlation")+ylab("transcriptomic correlation")+geom_text(data=as.data.frame(info), x=-0.5,y=0.5, aes(label = info))

}


plots = foreach(i = bin_md$V1) %dopar% {
  print(i)
  bin_plot(i)
}

plots_p = foreach(i = bin_md$V1) %dopar% {
  print(i)
  pbin_plot(i)
}

plots_c = foreach(i = bin_md$V1) %dopar% {
  print(i)
  cor_plot(i)
}


pdf("all_bins_pvals.pdf", onefile = TRUE)
for (i in seq(length(plots_p))) {
  print(i)
  do.call("grid.arrange", plots_p[i])
}
dev.off()

pdf("all_bins_FC.pdf", onefile = TRUE)
for (i in seq(length(plots))) {
  print(i)
  do.call("grid.arrange", plots[i])
}
dev.off()

pdf("all_bins_cor.pdf", onefile = TRUE)
for (i in seq(length(plots))) {
  print(i)
  do.call("grid.arrange", plots_c[i])
}
dev.off()


save(all_table,all_table_cor, dna_reads, rna_reads, cds_info, bin_md, station_md, file="main_var.Rdat")
