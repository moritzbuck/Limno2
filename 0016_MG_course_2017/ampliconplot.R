library(reshape2)
library(ggplot2)
library(plyr)
library(vegan)

# load the data

# load the OTU-table
table = read.table("OTU_table.txt", sep="\t", comment.char=";", head=T, row.names=1)

# load the predicted taxonomy with sintax
taxas = read.table("taxonomy.tax", sep="\t", as.is=T)
temp = taxas$V4
names(temp) = taxas$V1
taxas = temp

# extract the 2th taxonomical level, e.g. phylum
taxon.level = 2

taxas = sapply(strsplit(taxas,","),"[",taxon.level)

names(taxas) = sub(";size=[0-9]*;","", names(taxas))
# making a frequency table
freq_table = t(t(table)/colSums(table))

# manipulating to hell to get all in the right shape

mergy = merge(freq_table,as.data.frame(taxas), by="row.names")
mergy = mergy[,-1]
mergy = ddply(mergy, .(taxas), colwise(sum))
levels(mergy$taxas) = c(levels(mergy$taxas), "Unclassified", "Other")
cutoff = 0.01 # filtering taxas that are never more than `cutoff` percent from the data and calling them other
mergy$taxas[apply(mergy[,2:ncol(mergy)],1, max) < cutoff] = "Other"
mergy$taxas[is.na(mergy$taxas)] = "Unclassified"
mergy = melt(mergy, value.name="freq", variable.name="sample")


libraries = levels(as.factor(colnames(table)))

# extract metadata from library names
metadata=as.data.frame(t(as.data.frame(strsplit(sub("nasal_cavity","nasal",libraries),"_"), optional=FALSE)))

# make row.names and colnames
row.names(metadata) = libraries
colnames(metadata) = c("site","day")
metadata$day = as.numeric(metadata$day)

#  merging all in the enf
mergy=merge(mergy, metadata, by.x = 'sample', by.y="row.names")

# making "pretty" colors
ncols = length(levels(mergy$taxas))+1
cols = rainbow(ncols, s=.6, v=.9)[sample(1:ncols)]

#barplot
p1 = ggplot(mergy, aes(x=day, fill=taxas, y=freq))+geom_col()+facet_wrap(~site,ncol=1)+scale_fill_manual(values=cols)
ggsave("taxa_barplots.pdf", p1, width=12, height=6)

# compute Multidimensional scalings
dists.nasal = vegdist(t(table[,grepl("nasal",colnames(table))]))
dists.poop = vegdist(t(table[,!grepl("nasal",colnames(table))]))

mds.poop = data.frame(metaMDS(dists.poop)$points)
mds.poop = merge(mds.poop,metadata, by.x = 'row.names', by.y = 'row.names')
mds.nasal = data.frame(metaMDS(dists.nasal)$points)
mds.nasal = merge(mds.nasal,metadata, by.x = 'row.names', by.y = 'row.names')
# make pretty plots
p2 = ggplot(mds.nasal, aes(x=MDS1, y=MDS2, shape=site, col=day, label=Row.names))+scale_color_gradient(low="purple",high="green") + geom_label()
ggsave("MDS_nasal.pdf", p2)

p3 = ggplot(mds.poop, aes(x=MDS1, y=MDS2, shape=site, col=day, label=Row.names))+scale_color_gradient(low="purple",high="green") + geom_label()
ggsave("MDS_fecal.pdf", p3)
