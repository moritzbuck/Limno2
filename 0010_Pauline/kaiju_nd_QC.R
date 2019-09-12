library(plyr)
library(reshape2)
library(ggplot2)
library(vegan)
library(ggrepel)
library(dplyr)

barstak = function(level, proportion = TRUE, filter = TRUE, cutoff=0.01){
sub = ddply(kaiju_data[,c(samples,level)], c(level), colwise(sum))

rare = rowSums(sub[,samples])/sum(rowSums(sub[,samples])) < 0.01
sub[rare,level] = "Other"

if(filter) {
  sub=sub[sub[,level] != "Unclassified",]
}
if(proportion){
sub[,samples] = t(t(sub[,samples])/colSums(sub[, samples]))
}

mltd = melt(sub, value.name="abund", variable.name="library")

mltd = cbind(mltd, md[mltd$library,])
plot(ggplot(mltd, aes(x=as.factor(station), y=abund, fill=get(level)))+geom_col()+facet_grid(fract~type, scales="free_x"))
sub
}




# load the kaiju_data
mash = read.table("mash_distances.txt", as.is=TRUE)
md = read.csv("metadata.csv", h=TRUE, row.names=1, as.is =TRUE)
kaiju_data = read.csv("kaiju_table.csv",h=T, row.names=1, as.is = TRUE)
seq_data = read.table("seq_res.csv",h=T, sep=",")

seq_data$library = sub("_S[0-9_A-Z]*", "",seq_data$Sample)
seq_data$Dups = as.numeric(sub("%","",seq_data$Dups))
seq_data$GC = as.numeric(sub("%","",seq_data$GC))
seq_data$pair = sapply(strsplit(as.character(seq_data$Sample),"_"),"[",5)
seq_data$lane = sapply(strsplit(as.character(seq_data$Sample),"_"),"[",4)
seq_data$pair = sub("R([12])","\\1P", seq_data$pair, perl=T)
row.names(seq_data) = paste(seq_data$library, seq_data$pair, sep="_")
pre_dat = read.table("pre_filt_stats.tsv")
post_dat = read.table("post_filt_stats.tsv")
pre_dat$V1 = sub("R([12])","\\1P", pre_dat$V1, perl=T)
colnames(pre_dat) = c("library", "ori_volume")
colnames(post_dat) = c("library", "qced_volume")
volumes = merge(pre_dat, post_dat, by = "library")
volumes$kept = volumes$qced_volume/volumes$ori_volume
seq_data = merge(seq_data, volumes, by.x="row.names", by.y = "library")
row.names(seq_data) = seq_data$library.y
full.md = merge(seq_data, md, by.x="library", by.y="row.names")

#md = merge(md,seq_data, by.x = "row.names", by.y="library")


libraries = levels(as.factor(mash$V1))
libraries = sapply(strsplit(libraries,"/"), tail, 1)
libraries = tools::file_path_sans_ext(tools::file_path_sans_ext(libraries))
names(libraries) = levels(as.factor(mash$V1))
mash$V1 = libraries[mash$V1]
mash$V2 = libraries[mash$V2]


pmash = mash[,1:3]
distance_matrix = dcast(pmash, V1~V2)
row.names(distance_matrix) = distance_matrix$V1
distance_matrix = distance_matrix[,-1]

#pca = data.frame(predict(prcomp(distance_matrix))[,1:2])
pca = metaMDS(as.dist(distance_matrix))
pca = data.frame(pca$points)
pca$libraries = sub("_[12]P","",row.names(pca))
pca$pair = row.names(pca)
pca = merge(pca, md, by.x = "libraries", by.y="row.names")
pca = merge(pca, seq_data, by.x="pair", by.y="row.names")



samples = grep("P6404",colnames(kaiju_data), val=TRUE)
kaiju_data["0",is.na(kaiju_data["0",] )] = "Unclassified"
kaiju_data[is.na(kaiju_data)] = "UnclassifiedHere"


families = barstak("family", filter=FALSE, cutoff  = 0)
pca.taxo = data.frame(metaMDS(t(families[,-1]))$points)
pca.taxo = merge(pca.taxo, md, by.x = "row.names", by.y="row.names")


ggplot(pca.taxo, aes(x=MDS1, y=MDS2, label=libraries, shape=type, col=as.factor(fract), size=Dups)) + geom_label(size=0.1) + geom_point() + theme_bw()






maps = read.table("mapping_stats.csv", row.names=1)
colnames(maps)=c("rate")
maps$w_dups = read.table("mapping_stats_with_dups.csv", row.names=1)[row.names(maps),]
maps = merge(maps, md)
maps$libs = row.names(maps)

full.md$fullmap = read.table("full_coas_map.tsv", row.names=1)[full.md$library, "V2"]
full.md$binmap = read.table("good_bins_map.tsv", row.names=1)[full.md$library, "V2"]

full.md$fullmap[is.na(full.md$fullmap)] = 0
full.md$binmap[is.na(full.md$binmap)] = 0

full.md$eff.Mreads = full.md$Mreads*(1 - full.md$Dups/100)

data = read.table("bin_stats.csv", sep=",", head=T)
data$avg_prot_len = data$length*data$coding_density/data$nb_proteins
ggplot(data, aes(label=X, x=avg_prot_len, col=completeness, y=length))+geom_point()+scale_y_log10()+geom_text_repel(data = filter(data, avg_prot_len < 500 | avg_prot_len >1500))



nifH_data = read.table("nifHes.tsv", sep="\t", head=TRUE)
nifH_data = nifH_data[,!grepl(".var",colnames(nifH_data))]
melted.nifh = melt(nifH_data[,c(1,4:ncol(nifH_data))])
melted.nifh$type = sapply(melted.nifh$stations, function(x) md[md$station == x,"type"][1])
more_md = read.csv("../000_data/RNA_station_dat.csv", sep=",")
to_clean = names(which(apply(more_md, 2, function(x) sum(grepl("±", x))) > 0))
for(c in to_clean){
  more_md[,c] = sapply(strsplit(more_md[,c]," "), function(x) as.numeric(x[1]))
}
row.names(more_md) = sapply(strsplit(more_md$Station.nr," "),tail,1)
melted.nifh = cbind(melted.nifh,more_md[as.character(melted.nifh$stations),])
ggplot(melted.nifh, aes(x=stations, y=value, col=contigName))+geom_point()+geom_line()+facet_grid(~type)
