
library(ggplot2)
library(vegan)
library(reshape2)
library(rwantshue)
library(plyr)

input_tables = "../300_tables/unoise_sintax_silva/"
metadata_file = "../000_data/metadata.csv"
db_methyls = "../000_data/methylators.csv"

hgcA_tables = "../00A_hgcA_amplicon/otu_table_stump.csv"
hgcA_class = "../00A_hgcA_amplicon/classification_stump.txt"
sample_map = "../000_data/sample_map.csv"

table_file = file.path(input_tables,"otu_table.raw.csv")
family_file = file.path(input_tables,"otu_table.family.csv")

raw_otus = read.table(table_file, sep=",", ,h=T, comment.char="|", row.names=1)
raw_families = read.table(family_file, sep=",", ,h=T, comment.char="|", row.names=1)

ott = t(raw_otus[,1:(ncol(raw_otus)-1)])
curve = rarecurve(ott, step=300, label = FALSE)
dev.copy(pdf, "rarefaction_otus.pdf")
dev.off()

otf = t(raw_families[,1:(ncol(raw_families)-1)])
rarecurve(otf, step=300, label = FALSE)
dev.copy(pdf, "rarefaction_fams.pdf")
dev.off()
rared_otus = rarefy(ott[rowSums(ott) > 7,], 1692)
rared_families = rrarefy(otf[rowSums(otf) > 7,], 1692)

metadata = read.table(metadata_file,sep=",", h=T,row.names=1)
row.names(metadata) = paste("stubb", row.names(metadata), sep="")
metadata$GPS.N = as.numeric(as.vector(metadata$GPS.N))

hgca216s = read.table(sample_map, sep=",", ,h=F)
tt = hgca216s$V1
hgca216s = hgca216s$V2
hgca216s = paste("stubb", hgca216s, sep="")
names(hgca216s) = tt

raw_hgca_otus = read.table(hgcA_tables, sep=",", ,h=T, row.names=1)
row.names(raw_hgca_otus) = hgca216s[row.names(raw_hgca_otus)]

hgca_taxa = read.table(hgcA_class, sep="\t", row.names=1, as.is = T)

known_methyls = read.table(db_methyls, sep=",", h=T)
known_methyls$genus = sapply(strsplit(as.vector(known_methyls$Strain)," "),"[",1)

raw_taxas = raw_otus$taxa
names(raw_taxas) = row.names(raw_otus)

raw_otus = t(raw_otus[, 1:(ncol(raw_otus)-1)])

bad_samples = c("stubb200","stubb200")

filtered_otus = raw_otus[!rownames(raw_otus) %in% bad_samples,]
filtered_families = raw_families[!rownames(raw_families) %in% bad_samples,]
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

better_samples = row.names(metadata[!is.na(metadata$CS) & !is.na(metadata$CN) & !is.na(metadata$Water) & !is.na(metadata$Ratio.MeHg)  & row.names(metadata) %in% row.names(filtered_otus) ,])
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

family_shaped = melt(t(t(filtered_families)/colSums(filtered_families)))
family_shaped$Var1 = factor(family_shaped$Var1, levels=row.names(filtered_families)[order(rowSums(filtered_families), decreasing=T)])
family_shaped$pot_methyl = as.vector(family_shaped$Var1) %in%  levels(known_methyls$Family)
family_shaped$Name = metadata[family_shaped$Var2,"Name"]


scheme <- iwanthue(seed = 42, force_init = TRUE)
cols = scheme$hex(length(levels(family_shaped$Var1)), color_space = hcl_presets$pimp)

#hgca_taxa$bla = gsub(" \\([0-9]*\\)", "",sapply(strsplit(hgca_taxa$V2,";"),"[",4))
hgca_taxa$bla = gsub(" \\([0-9]*\\)", "",hgca_taxa$V2)

tt = as.data.frame(t(raw_hgca_otus/rowSums(raw_hgca_otus)))
tt$taxa = hgca_taxa[rownames(tt),"bla"]

hgca_shaped = melt(tt)
#hgca_shaped$Var1 = factor(hgca_shaped$Var1, levels=row.names(filtered_families)[order(rowSums(filtered_families), decreasing=T)])
hgca_shaped$Name = metadata[as.vector(hgca_shaped$variable),"Name"]

cols2 = scheme$hex(length(levels(as.factor(hgca_shaped$taxa))), color_space = hcl_presets$pimp)
1, levels=row.names(filtered_families)[order(rowSums(filtered_families), decreasing=T)])
family_shaped$pot_methyl = as.vector(family_shaped$Var1) %in%  levels(known_methyls$Family)
family_shaped$Name = metadata[family_shaped$Var2,"Name"]


scheme <- iwanthue(seed = 42, force_init = TRUE)
cols = scheme$hex(length(levels(family_shaped$Var1)), color_space = hcl_presets$pimp)

#hgca_taxa$bla = gsub(" \\([0-9]*\\)", "",sapply(strsplit(hgca_taxa$V2,";"),"[",4))
hgca_taxa$bla = gsub(" \\([0-9]*\\)", "",hgca_taxa$V2)

tt = as.data.frame(t(raw_hgca_otus/rowSums(raw_hgca_otus)))
tt$taxa = hgca_taxa[rownames(tt),"bla"]

hgca_shaped = melt(tt)
#hgca_shaped$Var1 = factor(hgca_shaped$Var1, levels=row.names(filtered_families)[order(rowSums(filtered_families), decreasing=T)])
hgca_shaped$Name = metadata[as.vector(hgca_shaped$variable),"Name"]

cols2 = scheme$hex(length(levels(as.factor(hgca_shaped$taxa))), color_space = hcl_presets$pimp)

hgca_taxa_2 = as.data.frame(t(sapply(strsplit(hgca_taxa$bla, ";"), function(x) {length(x)=7; return(x)})))
row.names(hgca_taxa_2) = row.names(hgca_taxa)
colnames(hgca_taxa_2) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
hgca_taxa_2$tip = NA
tt = apply(hgca_taxa_2,1, function(x) tail(x[!is.na(x)],1))
hgca_taxa_2[names(tt), "tip"] = as.vector(unlist(tt))

tt = as.data.frame(t(raw_hgca_otus) )

tt$taxa = hgca_taxa[row.names(tt),"bla"]
tt$taxa[tt$taxa== "" ] = NA
tt$taxa[is.na(tt$taxa) ] = "Unclassified"
write.csv(tt, "../300_tables/hgca_tables/hgca_table.raw.csv")

for(level in c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "tip"))
{
  tt = as.data.frame(t(raw_hgca_otus) )
  tt$taxa = hgca_taxa_2[row.names(tt),level]
  tt$taxa[tt$taxa== "" ] = NA
  levels(tt$taxa) = c(levels(tt$taxa), "Unclassified")
  tt$taxa[is.na(tt$taxa) ] = "Unclassified"
  tt = ddply(tt, .(taxa), colwise(sum))
  write.csv(tt, sprintf("../300_tables/hgca_tables/hgca_table.%s.csv", level))
}


rared_families = rared_families[!row.names(rared_families) %in% bad_samples,]
öb = row.names(metadata[metadata$Name == "Strömsjöliden",])
öb = öb [öb %in% row.names(rared_families)]

mds=metaMDS(rared_families[öb,], try = 100, trymax=100)
mds_points = mds$points
mds_points = cbind(mds_points, metadata[row.names(mds_points),])

bla = envfit(mds, metadata[row.names(mds_points),c("MeHg","Ratio.MeHg","Thg","Water","CN","CS", "Name")], na.rm = TRUE)
arrows = as.data.frame(bla$vector$arrow)
arrows = arrows*t(bla$vector$r)

arrows2 = as.data.frame(bla$factors$centroids)

arrows = rbind(arrows, arrows2)

arrows$label = row.names(arrows)

ggplot(mds_points, aes(MDS1, MDS2, col=Name, size=Ratio.MeHg))+geom_point()+geom_segment(data=arrows, aes(xend=NMDS1,yend=NMDS2), x=0, y=0, col="red", size=1) + geom_text(data=arrows, aes(x=NMDS1,y=NMDS2, label =label), col="black", size=8)


anaeros = cbind(reads = rared_families[,"Anaerolineaceae"], metadata[row.names(rared_families),])
ggplot(anaeros, aes(x=Ratio.MeHg, y=reads/1692, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Anaerolineaceae")
dev.copy(pdf,"anaeros.pdf", width=8, height=2.5)
dev.off()

anaeros = cbind(reads = rared_families[,"Syntrophobacteraceae"], metadata[row.names(rared_families),])
ggplot(anaeros, aes(x=Ratio.MeHg, y=reads/1692, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Syntrophobacteraceae")
dev.copy(pdf,"syntrophos.pdf", width=8, height=2.5)
dev.off()

anaeros = cbind(reads = rared_families[,"Desulfobulbaceae"], metadata[row.names(rared_families),])
ggplot(anaeros, aes(x=Ratio.MeHg, y=reads/1692, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Desulfobulbaceae")
dev.copy(pdf,"desulfos.pdf", width=8, height=2.5)
dev.off()

anaeros = cbind(reads = rared_families[,"Geobacteraceae"], metadata[row.names(rared_families),])
ggplot(anaeros, aes(x=Ratio.MeHg, y=reads/1692, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Geobacteraceae")
dev.copy(pdf,"geos.pdf", width=8, height=2.5)
dev.off()



level = "Family"

  tt = as.data.frame(t(raw_hgca_otus) )
  tt$taxa = hgca_taxa_2[row.names(tt),level]
  tt$taxa[tt$taxa== "" ] = NA
  levels(tt$taxa) = c(levels(tt$taxa), "Unclassified")
  tt$taxa[is.na(tt$taxa) ] = "Unclassified"
  tt = ddply(tt, .(taxa), colwise(sum))

 hgca_mds = metaMDS(t(tt[,2:ncol(tt)]), try = 100, trymax=100)
 hgca_mds_points = hgca_mds$points
 hgca_mds_points = cbind(hgca_mds_points, metadata[row.names(hgca_mds_points),])

 bla = envfit(hgca_mds, metadata[row.names(hgca_mds_points),c("MeHg","Ratio.MeHg","Thg","Water","CN","CS", "Name")], na.rm = TRUE)
 arrows = as.data.frame(bla$vector$arrow)
 arrows = arrows*t(bla$vector$r)

 arrows2 = as.data.frame(bla$factors$centroids)

 arrows = rbind(arrows, arrows2)

 arrows$label = row.names(arrows)

 ggplot(hgca_mds_points, aes(MDS1, MDS2, col=Name, size=Ratio.MeHg))+geom_point()+geom_segment(data=arrows, aes(xend=NMDS1,yend=NMDS2), x=0, y=0, col="red", size=1) + geom_text(data=arrows, aes(x=NMDS1,y=NMDS2, label =label), col="black", size=8)

tt2 = t(t(tt2)/colSums(tt2))
anaeros = cbind(reads = tt2[tt$taxa == "Ruminococcaceae",], metadata[colnames(tt2),])

ggplot(anaeros, aes(x=CS, y=reads, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Ruminococcaceae")
dev.copy(pdf,"ruminos_hgca.pdf", width=8, height=2.5)
dev.off()

anaeros = cbind(reads = tt2[tt$taxa == "Desulfuromonadaceae",], metadata[colnames(tt2),])

ggplot(anaeros, aes(x=CS, y=reads, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Desulfuromonadaceae")
dev.copy(pdf,"desulfos_hgca.pdf", width=8, height=2.5)
dev.off()

anaeros = cbind(reads = tt2[tt$taxa == "Methanomassiliicoccaceae",], metadata[colnames(tt2),])

ggplot(anaeros, aes(x=CS, y=reads, col=Name))+geom_point()+facet_wrap(~Name, scale="free_x")+ylab("relative abundance")+ggtitle("Methanomassiliicoccaceae")
dev.copy(pdf,"methanomass_hgca.pdf", width=8, height=2.5)
dev.off()
