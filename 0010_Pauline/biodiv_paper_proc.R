library(data.table)
library(ggplot2)
library(qualpalr)
library(vegan)
library(ggrepel)

raw_reads = fread("read_counts_all_samples.csv")

md = read.csv("metadata.csv")
sample_map = read.csv("sample_map.csv")
sample_map = merge(sample_map, read.table("lib_sizes.txt"), by.x="X", by.y="V1")

library_md = merge(sample_map, md, by.x="sample", by.y="X")
lib_sizes = merge(melt(raw_reads), library_md, by.x = "variable", by.y = "X")[, .(lib_sizes = sum(value)), by = .(sample)]
setkey(lib_sizes, "sample")
data = melt(raw_reads, id.vars="V1", variable.name = "library", value.name = "raw_reads")
colnames(data)[1] = "MAB"


mab_md = fread("all_good_mags.stats")
mab_md$type = "coass"
mab_md$type[grepl("P6404",mab_md$bin_name) |  grepl("Sample",mab_md$bin_name)] = "single"

mab_md[,mab_id := paste(type, assmbler, bin_name, sep="_")]
mab_md[,mab_id := gsub("-","_",mab_id)]

setkey(mab_md, mab_id)

tax = fread("sourmash_lca_gtdbtk_combo_plus_dirty_gtdbtk.tax", sep=",")
tax[,MAGs := gsub("-","_",MAGs)]

mab_md = merge(mab_md, tax, by.x = "mab_id", by.y = "MAGs", all.x = TRUE )

viral_dat = fread("viral_fraction.csv"); setkey(viral_dat, mab_id)
mab_md = merge(mab_md, viral_dat, by = "mab_id", all.x = TRUE )
mab_md[, viral_fract := viral_length/length]

data[, MAB := gsub("-","_",MAB)]

data = merge(data, library_md, by.x = "library", by.y = "X")
data$read_fract = data$raw_reads/lib_sizes[as.character(data$sample)]$lib_sizes

data = merge(data, mab_md, by.x = "MAB", by.y = "mab_id", all.x = TRUE)
bla = data[,.(read_fract=mean(read_fract), max_fract = max(read_fract), brine_fract = mean(read_fract[type.x == "Brine"]), water_fract = mean(read_fract[type.x == "Seawater"]) ) , by = .(MAB)]
setkey(bla,MAB)
mab_md[,tot_cov := bla[mab_md]$read_fract]
mab_md[, brine_cov :=  bla[mab_md]$brine_fract]
mab_md[, Seawater_cov :=  bla[mab_md]$water_fract]
mab_md[,max_cov := bla[mab_md]$max_fract]

data$set = "MAB"
data[!is.na(viral_fract) & viral_fract > 0.1, Phylum := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, Class := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, Order := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, Family := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, set := "viral" ]
data[MAB == "*", set := "unassembled"]
data[MAB == "unbinned", set := "unbinned"]



rare_proc = 0.01
temp = data[, .(read_fract = sum(read_fract), fract = fract[1], station = station[1], type = type.x[1]), by = .(Domain, Phylum, Class, Order, Family,Genus, library, set)]

temp = temp[, tax_str := paste(Domain,Phylum, Class,Order,Family, sep = ";" ) ]
temp[set == 'viral', tax_str := "viral"]
temp = temp[set != 'unbinned' & set != 'unassembled' & tax_str != "NA;NA;NA;NA;NA"]
temp[, read_fract := read_fract/sum(read_fract) ,by = .(station, type, fract)]
rares = temp[, .(read_fract = sum(read_fract)),by = .(tax_str, station, type, fract)][, .(max_fract = max(read_fract)), by = tax_str][max_fract < rare_proc]$tax_str
temp[tax_str %in% rares, tax_str := "rares"]

tt = sapply(strsplit(levels(factor(temp$tax_str)),";"), function(x) paste(head(x,3),sep=";", collapse=";"))

tax_groups = sapply(levels(factor(tt)), function(x) grep(x, levels(factor(temp$tax_str)), val=TRUE) )
tax_groups = tax_groups[sapply(tax_groups, length) > 0 ]
start = qualpal(length(tax_groups)-5, list(h=c(0,360), s = c(0.7,0.7), l=c(0.6, 0.6)))$HSL[,1]
names(start) = names(tax_groups)[1:(length(tax_groups)-5)]

cols = lapply(names(tax_groups)[1:(length(tax_groups)-5)], function(x) if(length(tax_groups[[x]]) > 1) qualpal(n = length(tax_groups[[x]]), list(h = c(start[x], start[x]), s = c(0.6, 0.9), l = c(0.6, 0.85)))$hex else hsv(start[[x]]/360, 0.8,0.8) )

greys = qualpal(n=3, list(h=c(200,200), s=c(0,0), l=c(0.4, 0.85)))$hex

cols = c(unlist(cols), rev(greys))

#tax_groups = sapply(levels(factor(paste(temp$Domain,temp$Phylum,temp$Class, sep=";"))), function(x) grep(x, levels(factor(temp$tax_str)), val=TRUE) )
#tax_groups = tax_groups[sapply(tax_groups, length) > 0 ]
#start = qualpal(length(tax_groups), list(h=c(0,360), s = c(0.7,0.7), l=c(0.6, 0.6)))$HSL[,1]
#names(start) = names(tax_groups)
#cols = lapply(names(start), function(x) if(length(tax_groups[[x]]) > 1) qualpal(n = length(tax_groups[[x]]), list(h = c(start[x], start[x]), s = c(0.6, 0.9), l = c(0.6, 0.85)))$hex else hsv(start[[x]]/360, 0.7,0.6) )
#greys = qualpal(n=2, list(h=c(200,200), s=c(0,0), l=c(0.4, 0.85)))$hex

#cols = c(unlist(cols), greys)

fam_data = temp


families_barstack = ggplot(temp, aes(x = factor(station), y=read_fract, fill = tax_str))+theme_minimal()+ geom_bar(stat = 'identity')+facet_wrap(type~factor(fract), scales = "free")+guides(fill=guide_legend(ncol=2))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.text=element_text(size=7), legend.key.size = unit(0.5,"line"))+guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=cols)



temp = data[, .(read_fract = sum(read_fract), fract = fract[1], station = station[1], type = type.x[1]), by = .(Domain, Phylum, Class, Order, Family,Genus, library, set)]

temp = temp[, tax_str := paste(Domain,Phylum, Class, sep = ";" ) ]
temp$tax_str[temp$tax_str == "NA;NA;NA"] = "unclassified MAB"
temp$tax_str[temp$set == 'unbinned']  = "unbinned"
temp$tax_str[temp$set == 'unassembled']  = "unassembled"
temp$tax_str[temp$set == 'viral']  = "viral"
rares = temp[, .(read_fract = sum(read_fract)),by = .(tax_str, station, type, fract)][, .(max_fract = max(read_fract)), by = tax_str][max_fract < rare_proc]$tax_str
temp[tax_str %in% rares, tax_str := "rares"]
temp[tax_str %in% rares, Phylum := NA]
temp[tax_str %in% rares, Class := NA]
temp[tax_str %in% rares, Domain := NA]

tt = sapply(strsplit(levels(factor(temp$tax_str)),";"), function(x) paste(head(x,2),sep=";", collapse=";"))

tax_groups = sapply(levels(factor(tt)), function(x) grep(x, levels(factor(temp$tax_str)), val=TRUE) )
tax_groups = tax_groups[sapply(tax_groups, length) > 0 ]
start = qualpal(length(tax_groups)-5, list(h=c(0,360), s = c(0.7,0.7), l=c(0.6, 0.6)))$HSL[,1]
names(start) = names(tax_groups)[1:(length(tax_groups)-5)]

cols = lapply(names(tax_groups)[1:(length(tax_groups)-5)], function(x) if(length(tax_groups[[x]]) > 1) qualpal(n = length(tax_groups[[x]]), list(h = c(start[x], start[x]), s = c(0.6, 0.9), l = c(0.6, 0.85)))$hex else hsv(start[[x]]/360, 0.8,0.8) )

greys = qualpal(n=5, list(h=c(200,200), s=c(0,0), l=c(0.4, 0.85)))$hex

cols = c(unlist(cols), rev(greys))

coarse_barstack = ggplot(temp, aes(x = factor(station), y=read_fract, fill = tax_str))+ geom_bar(stat = 'identity')+theme_minimal()+facet_wrap(type~factor(fract), scales = "free")+guides(fill=guide_legend(ncol=2))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.text=element_text(size=7), legend.key.size = unit(0.5,"line"))+guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=cols)

pfam_matrix = read.csv("pfam_matrix.csv", row.names=1, )
pfam_matrix = pfam_matrix[rowSums(pfam_matrix) != 0,]

full.pfam.mds = metaMDS(t(pfam_matrix))
pfam.species = as.data.frame(full.pfam.mds$species)
pfam.mds = as.data.frame(full.pfam.mds$point)
pfam.mds$lib = gsub(".","-",sub(".bam","",row.names(pfam.mds)), fixed=TRUE)

pfam.mds = merge(merge(pfam.mds, sample_map, by.x = "lib", by.y = "X"), md, by.x = "sample", by.y="X")
bigs = pfam.species[sqrt(pfam.species$MDS1**2 + pfam.species$MDS2**2) > 1 & !is.na(pfam.species$MDS1),]
bigs$name = row.names(bigs)

pfam_mds = ggplot()+geom_point(data = pfam.species, aes(x=MDS1, y=MDS2), col="purple", shape=20, alpha=0.1)+ geom_label_repel(data = bigs, aes(x=MDS1, y=MDS2, label = name)) + geom_point(data = pfam.mds, aes(x=MDS1, y=MDS2, col = as.factor(fract), shape=type),size = 5, alpha = 1)


weird = c("P6404_212","P6404_234","P6404_267","P6404_282","P6404_283","P6404_285","P6404_289","P6404_290","P6404_291")

freq_mat_families = dcast(data = fam_data[!library %in% weird & type == "Brine"], formula = tax_str~library, value.var = "read_fract", fun.aggregate=sum)

full.tax.mds = metaMDS(t(freq_mat_families[,-1]))
tax.species = as.data.frame(full.tax.mds$species)
#tax.species$name = sapply(strsplit(freq_mat_families$tax_str,";"), tail,1)
tax.species$name = freq_mat_families$tax_str

tax.mds = as.data.frame(full.tax.mds$point)
tax.mds$lib = colnames(freq_mat_families[,-1])

tax.mds = merge(merge(tax.mds, sample_map, by.x = "lib", by.y = "X"), md, by.x = "sample", by.y="X")
bigs = tax.species[sqrt(tax.species$MDS1**2 + tax.species$MDS2**2) > 0.7,]

brine_mds = ggplot()+geom_point(data = tax.species, aes(x=MDS1, y=MDS2), col="purple", shape=20, alpha=0.5, size=3)+ geom_label_repel(data = bigs, aes(x=MDS1, y=MDS2, label = name)) + geom_point(data = tax.mds, aes(x=MDS1, y=MDS2, col = as.factor(fract), shape=type),size = 5, alpha = 1)

freq_mat_families = dcast(data = fam_data[!library %in% weird & type == "Seawater"], formula = tax_str~library, value.var = "read_fract", fun.aggregate=sum)


full.tax.mds = metaMDS(t(freq_mat_families[,-1]))
tax.species = as.data.frame(full.tax.mds$species)
#tax.species$name = sapply(strsplit(freq_mat_families$tax_str,";"), tail,1)
tax.species$name = freq_mat_families$tax_str

tax.mds = as.data.frame(full.tax.mds$point)
tax.mds$lib = colnames(freq_mat_families[,-1])

tax.mds = merge(merge(tax.mds, sample_map, by.x = "lib", by.y = "X"), md, by.x = "sample", by.y="X")
bigs = tax.species[sqrt(tax.species$MDS1**2 + tax.species$MDS2**2) > 0.7,]

sea_mds = ggplot()+geom_point(data = tax.species, aes(x=MDS1, y=MDS2), col="purple", shape=20, alpha=0.5, size=3)+ geom_label_repel(data = bigs, aes(x=MDS1, y=MDS2, label = name)) + geom_point(data = tax.mds, aes(x=MDS1, y=MDS2, col = as.factor(fract), shape=type),size = 5, alpha = 1)


nif = read.csv("nifCovs.csv", h=F, col.names=c("lib","cov"))
nif$lib = sub(".","-", sub(".bam","",nif$lib), fixed=TRUE)
nif = merge(merge(nif, sample_map, by.x = "lib", by.y = "X"), md, by.x = "sample", by.y="X")

cyano_ags = read.table("cyanophycin_mags.txt", as.is=TRUE)$V1
logic =!grepl("Sample", cyano_ags) & !grepl("P64", cyano_ags)
cyano_ags[logic] =  sub("single_", "coass_", cyano_ags[logic])

photo = read.table("photosystem_mags.txt", as.is = TRUE)
photo = photo$V2[photo$V1 > 15]
logic =!grepl("Sample", photo) & !grepl("P64", photo)
photo[logic] =  sub("single_", "coass_", photo[logic])

plastocyanin = read.table("plastocyanin_mags.txt", as.is = TRUE)
plastocyanin = plastocyanin$V2[plastocyanin$V1 > 0]
logic =!grepl("Sample", plastocyanin) & !grepl("P64", plastocyanin)
plastocyanin[logic] =  sub("single_", "coass_", plastocyanin[logic])

pot_euks  = c('coass_megahit_Station_3_246',
 'single_megahit_P6404_225_178',
 'coass_megahit_FullCoassembly_4958',
 'coass_megahit_Station_2_969',
 'coass_megahit_big_fraction_2137',
 'coass_megahit_FullCoassembly_4622',
 'coass_megahit_FullCoassembly_3841',
 'coass_megahit_big_fraction_3805',
 'coass_megahit_FullCoassemblyNo3_1901',
 'coass_megahit_Station_1_162',
 'coass_megahit_Station_9_200',
 'coass_megahit_FullCoassemblyNo3_2449',
 'coass_megahit_FullCoassemblyNo3_3545',
 'coass_megahit_FullCoassemblyNo3_553',
 'coass_megahit_big_fraction_1890',
 'coass_megahit_Station_9_398',
 'coass_megahit_FullCoassemblyNo3_2721',
 'coass_megahit_Station_3_438',
 'single_megahit_P6404_249_63',
 'coass_megahit_FullCoassembly_4732',
 'coass_megahit_FullCoassembly_2236',
 'coass_megahit_FullCoassembly_4348')

pot_euk_tax = c('Eukaryota(1.00);Stramenopiles(0.98);Bacillariophyta(0.98);Bacillariophyceae(0.66);Bacillariophycidae(0.66)', 'Eukaryota(1.00);Alveolata(0.71);Ciliophora(0.71);Intramacronucleata(0.71);Spirotrichea(0.71);Stichotrichia(0.71);Sporadotrichida(0.71);Oxytrichidae(0.71)', 'Eukaryota(0.89);Stramenopiles(0.61)', 'Eukaryota(1.00);Alveolata(0.83);Ciliophora(0.83);Intramacronucleata(0.83);Spirotrichea(0.83);Stichotrichia(0.83);Sporadotrichida(0.83);Oxytrichidae(0.83)', 'Eukaryota(0.82)', 'Eukaryota(0.82)', 'Eukaryota(1.00);Alveolata(0.71);Ciliophora(0.71);Intramacronucleata(0.71);Spirotrichea(0.71);Stichotrichia(0.71);Sporadotrichida(0.71);Oxytrichidae(0.71)', 'Eukaryota(1.00);Alveolata(0.71);Ciliophora(0.71);Intramacronucleata(0.71);Spirotrichea(0.71);Stichotrichia(0.71);Sporadotrichida(0.71);Oxytrichidae(0.71)', 'Eukaryota(0.85);Stramenopiles(0.50)', 'Eukaryota(0.94);Viridiplantae(0.94);Chlorophyta(0.88);Mamiellophyceae(0.88);Mamiellales(0.88);Mamiellaceae(0.88);Micromonas(0.88);Micromonas commoda(0.59)', 'Eukaryota(0.75);Stramenopiles(0.42)', 'Eukaryota(0.81);Opisthokonta(0.44)', 'Eukaryota(0.82)', 'Eukaryota(0.75);Haptophyceae(0.54)', 'Eukaryota(0.98);Stramenopiles(0.96);Bacillariophyta(0.96);Bacillariophyceae(0.62);Bacillariophycidae(0.62)', 'Eukaryota(0.93);Stramenopiles(0.79)', 'Eukaryota(0.84);Stramenopiles(0.48)', 'Eukaryota(0.88)', 'Eukaryota(1.00);Stramenopiles(0.98);Bacillariophyta(0.98);Bacillariophyceae(0.65);Bacillariophycidae(0.65)', 'Eukaryota(0.82)', 'Eukaryota(0.79);Stramenopiles(0.64)', 'Eukaryota(0.98);Stramenopiles(0.96);Bacillariophyta(0.96);Bacillariophyceae(0.64);Bacillariophycidae(0.64)')

buscos = read.table("busco_table.tsv", row.names=1)


names(pot_euk_tax) = pot_euks
euks_md = mab_md[pot_euks]
euks_md[, c("viral_fract","viral_length", "Domain", "Phylum", "Class", "Order", "Family", "Genus","Species", "bin_name", "type", "assmbler", "file", "strain_heterogeneity", "taxo:checkm", "V1") := NULL]
euks_md[, "putative_tax" := pot_euk_tax[mab_id]]
euks_md[, buscos_completness := buscos[mab_id,]]

pps = c(coarse_barstack, families_barstack, brine_mds, sea_mds, pfam_mds)
pps = c(coarse_barstack, families_barstack, brine_mds, sea_mds, pfam_mds)
ggsave("coarse_barstack.pdf",coarse_barstack, width=15, height=12)
ggsave("families_barstack.pdf",families_barstack, width=15, height=12)
ggsave("brine_mds.pdf",brine_mds, width=15, height=12)
ggsave("sea_mds.pdf",sea_mds, width=15, height=12)
ggsave("pfam_mds.pdf",pfam_mds, width=15, height=12)
