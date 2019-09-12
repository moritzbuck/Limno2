library(ggplot2)
library(vegan)

input_tables = "../300_tables/unoise_sintax_silva/"
metadata_file = "../000_data/metadata.csv"
db_methyls = "../000_data/methylators.csv"

hgcA_tables = "../00A_hgcA_amplicon/otu_table_stump.csv"
sample_map = "../000_data/sample_map.csv"
table_file = file.path(input_tables,"otu_table.raw.csv")

raw_otus = read.table(table_file, sep=",", ,h=T, comment.char="|", row.names=1)
metadata = read.table(metadata_file,sep=",", h=T,row.names=1)
row.names(metadata) = paste("stubb", row.names(metadata), sep="")

hgca216s = read.table(sample_map, sep=",", ,h=F)
tt = hgca216s$V1
hgca216s = hgca216s$V2
hgca216s = paste("stubb", hgca216s, sep="")
names(hgca216s) = tt

raw_hgca_otus = read.table(hgcA_tables, sep=",", ,h=T, row.names=1)

known_methyls = read.table(db_methyls, sep=",", h=T)
known_methyls$genus = sapply(strsplit(as.vector(known_methyls$Strain)," "),"[",1)
sapply(strsplit(as.vector(raw_taxas),";"),"[", 3)


raw_taxas = raw_otus$taxa
names(raw_taxas) = row.names(raw_otus)

raw_otus = t(raw_otus[, 1:(ncol(raw_otus)-1)])

bad_samples = c("stubb200","stubb200")

filtered_otus = raw_otus[!rownames(raw_otus) %in% bad_samples,]
filtered_taxas = raw_taxas[colnames(filtered_otus)]

fams = sapply(strsplit(as.vector(filtered_taxas),";"),"[", 5)
pot_methyl = names(filtered_taxas[which(fams %in% levels(known_methyls$Family))])

mds=metaMDS(filtered_otus, try = 100, trymax=100)

mds_points = mds$points
mds_points = cbind(mds_points, metadata[row.names(mds_points),])


pot_filtered_otus = filtered_otus[,pot_methyl]
pot_filtered_otus = pot_filtered_otus[rowSums(pot_filtered_otus) > 10,]
pot_mds=metaMDS(pot_filtered_otus, try = 100, trymax=100)

pot_mds_points = pot_mds$points
pot_mds_points = cbind(pot_mds_points, metadata[row.names(pot_mds_points),])

better_samples = row.names(metadata[!is.na(metadata$CS) & !is.na(metadata$CN) & !is.na(metadata$Water) & !is.na(metadata$Ratio.MeHg) & row.names(metadata) %in% row.names(filtered_otus),])
more_filtered_otus = filtered_otus[better_samples,]

adonis2(more_filtered_otus~Name+Water+Thg+Ratio.MeHg+CN+CS+treat_correct, metadata[row.names(more_filtered_otus),])

hgca_mds=metaMDS(raw_hgca_otus, try = 100, trymax=100)

hgca_mds_points = hgca_mds$points
hgca_mds_points = cbind(hgca_mds_points, metadata[row.names(hgca_mds_points),])

sub_mds=metaMDS(filtered_otus[row.names(raw_hgca_otus),], try = 100, trymax=100)

sub_mds_points = sub_mds$points
sub_mds_points = cbind(sub_mds_points, metadata[row.names(sub_mds_points),])

mantel(vegdist(raw_hgca_otus), vegdist(filtered_otus[row.names(raw_hgca_otus),]))

methyl_frac=cbind(fraction = rowSums(filtered_otus[,pot_methyl]) / rowSums(filtered_otus) , metadata[row.names(filtered_otus),])
