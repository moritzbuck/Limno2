library(data.table)
library(dada2)
seqs = getSequences("all_matam.fasta")
taxo = assignTaxonomy(seqs, refFasta="~/dbs/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxo = addSpecies(taxo, "~/dbs/silva_species_assignment_v132.fa.gz")
taxo = as.data.table(taxo)
taxo$id = sapply(strsplit(names(seqs)," "),head,1)
taxo$count = sapply(strsplit(names(seqs)," "),function(x) as.numeric(sub("count=", "",x[3])))
taxo$sample = sapply(strsplit(taxo$id,"_"), function(x) paste(x[1:(length(x)-1)], collapse="_"))
setkey(taxo,"id")
md = read.csv("9999_metadata/metadata.csv")
sample_map = read.csv("9999_metadata/sample_map.csv", row.names = 1)
taxo$biolo = sample_map[taxo$sample,"sample"]
taxo = merge(taxo, md, by.x = "biolo", by.y = "X" )
write.table(taxo[,.(id, tax_str = apply(.SD[,.(Kingdom, Phylum, Class, Order, Family, Genus, Species)],1, function(x) paste(x[!is.na(x)], sep=";", collapse=";")))],"matam_mapping.tax", row.names=FALSE, col.names=FALSE, quote=FALSE)
lib_sizes = fread("1000_processed_reads/lib_sizes.txt")
setkey(lib_sizes, V1)
lib_sizes$sample = sample_map[lib_sizes$V1,"sample"]
biolo_counts = lib_sizes[, .(reads = sum(V2)), by = sample]
setkey(biolo_counts, sample)
taxo$biolo_size = biolo_counts[as.character(taxo$biolo), reads]
taxo$norm_reads = taxo$count/taxo$biolo_size
freqs = dcast(taxo [, .( freq = sum(norm_reads) ), by = .(paste(Kingdom,Phylum,Class,Order, sep=";") , biolo)], paste~biolo, fill=0)
euks = as.data.frame(freqs[grepl("Eukar", paste)])
row.names(euks) = euks[,1]
euks = euks[,-1]
write.table(euks, "euks_table.csv")
write.table(taxo, "full_matam.table.csv")

mds = metaMDS(t(euks[,intersect(as.character(md[md$fract == 3.0, "X"]), colnames(euks))]), try=100)
spec = mds$species
spec = as.data.frame(spec[spec[,1] == spec[,1] & ! is.na(spec[,1]) ,])
spec$phylum = sapply(strsplit(row.names(spec),";"),"[", 2)

cols = c("#cefff4","#e224df","#3c7a00","#c954ff","#01ac64","#ff77e8","#ffbd3e","#001241","#ffd88d","#84007f","#d5ffd3","#ff2b3a","#0195ca","#d36800","#9892ff","#a91e00","#b4ceff","#122000","#ffa1cc","#001312","#a5005d","#005c4f","#014670")

ggplot(data = spec , aes(x=MDS1, y = MDS2, col=phylum)) + geom_point(shape=16, size=2) + geom_point( data = merge(mds$points,md, by.x="row.names", by.y = "X"), aes(shape=type), col="black", size=4)+ geom_label_repel( data = merge(mds$points,md, by.x="row.names", by.y = "X"), aes(label=Row.names), col="black", size=2)+ scale_shape_manual(values=c(17,15))+ scale_color_manual(values=cols)+theme_bw()

ggsave("large_fract_euks_mds.pdf", width=12, height=8)

mds = metaMDS(t(euks[,intersect(as.character(md[md$fract == 0.8 & md$X != "P6404_210", "X"]), colnames(euks))]), try=100)
spec = mds$species
spec = as.data.frame(spec[spec[,1] == spec[,1] & ! is.na(spec[,1]) ,])
spec$phylum = sapply(strsplit(row.names(spec),";"),"[", 2)

cols = c("#da4e8d","#61b650","#a44ec9","#b1b23a","#6164d3","#de9035","#b87dde","#627c35","#d052b2","#4cb08d","#d33d55","#4aadd6","#d4532e","#6782ca","#be9b59","#7c509b","#9f5d2e","#d394d1","#a1475a","#e0807f","#9e5284")

ggplot(data = spec , aes(x=MDS1, y = MDS2, col=phylum)) + geom_point(shape=16, size=2) + geom_point( data = merge(mds$points,md, by.x="row.names", by.y = "X"), aes(shape=type), col="black", size=4)+ geom_label_repel( data = merge(mds$points,md, by.x="row.names", by.y = "X"), aes(label=Row.names), col="black", size=2)+ scale_shape_manual(values=c(17,15))+ scale_color_manual(values=cols)+theme_bw()

ggsave("medium_fract_euks_mds.pdf", width=12, height=8)
