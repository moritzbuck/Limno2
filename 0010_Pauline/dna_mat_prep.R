library(data.table)
library(vegan)
library(ggplot2)

map_raw = fread("2000_MAG_sets/all_mags/cat/mapping/map_table.tsv")
map_raw = map_raw[,grep("-var",colnames(map_raw), invert=TRUE, value=TRUE),with = FALSE]
melted_raw = melt(map_raw,  measure.vars=grep("bam",colnames(map_raw),value=TRUE))
mag_stats = fread("2000_MAG_sets/all_mags/all_good_mags.stats")
setkey(mag_stats, V1)
tt = gsub("-","_",map_raw$contigName)
tt = strsplit(tt,"_")
tt = sapply(tt, function(x) paste(x[2:(length(x)-1)], collapse="_"))
names(tt) = map_raw$contigName
melted_raw[,bin := tt[contigName]]
binWiseCov = melted_raw[, .(cov = sum(contigLen*value)), by = .(bin,variable)]

binWiseCov$variable = sub(".bam","",binWiseCov$variable)
binWiseCov = merge(binWiseCov, sample_map, by.x="variable", by.y="V1")
binWiseReads = binWiseCov[ , .(reads = sum(cov)) , by = .(V2, bin)]


lib_sizes = fread("1000_processed_reads/lib_sizes.txt")
setkey(lib_sizes, V1)
lib_sizes[,coef := V2/mean(V2)]
map_rates = fread("2000_MAG_sets/all_mags/cat/mapping/mapping_rates.txt",header=FALSE, key="V1")

good_reads = as.data.table(read.table("2000_MAG_sets/all_mags/cat_goods/mapping/mapping_rates.txt", h=TRUE))
setkey(good_reads, mated)

map_rates = merge(map_rates, good_reads, by.x="V1", by.y="sample")
map_rates[, unmapped := 100-V2]
map_rates[,fwd := NULL]
map_rates[,rev := NULL]
map_rates[, unbinned := V2-mated]


taxo = fread("2000_MAG_sets/all_mags/good_genomes.gtdbtk.tax", fill = TRUE, col.names=c("bin","domain","phylum", "class","order","family", "genus", "species"))
taxo$bin = gsub("single_","",taxo$bin)
setkey(taxo, bin)

binWiseReads = merge(binWiseReads, taxo, by = "bin", all.x=TRUE)

md = fread("9999_metadata/metadata.csv")
sample_map = fread("9999_metadata/sample_map.csv")
md = merge(sample_map, md, by.x = "V2", by.y = "V1")

merge_table = function(level){
familyWise = binWiseReads[, .(reads = sum(reads)), by = .(V2,paste(domain,phylum, class, order, family, sep=";"))]

castFamilyWise = as.data.frame(dcast(familyWise, paste~V2))

row.names(castFamilyWise) = castFamilyWise$paste
castFamilyWise = castFamilyWise[,-1]

tot_covs = colSums(castFamilyWise)
tt=map_rates[names(tot_covs),V2]
tt[is.na(tt)] = mean( tt, na.rm=TRUE)

castFamilyWise["unassembled",] = tot_covs*tt[,unmapped]


castFamilyWise = t(t(castFamilyWise)/colSums(castFamilyWise))
familyWise = as.data.table(melt(castFamilyWise, value.name="coverage", c("taxon", "library")))

familyWise = merge(familyWise, md[!duplicated(V2)],by.x="library", by.y="V2")

mds = metaMDS(t(castFamilyWise))

points = merge(mds$points, md, by.x ="row.names", by.y = "V2")
species = data.frame(mds$species)
species$taxon = row.names(species)


ggplot(species, aes(x=MDS1, y=MDS2, label=taxon))+geom_text(size=1)+geom_point(data = points, aes(col=factor(fract), label=NULL, shape=type))


ggsave("mds_genomic.pdf")
}
