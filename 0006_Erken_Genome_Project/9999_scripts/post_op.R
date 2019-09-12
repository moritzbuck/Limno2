library(vegan)

data = read.table("../3000_tables/otu_table.raw.csv",sep = ",", as.is=TRUE, header=TRUE, comment.char="$", row.names=1)
counts = data[,-ncol(data)]
sample_map = read.csv("../sample_map.csv", sep="\t", row.names = 1, as.is = TRUE)
row.names(sample_map) = sub("P1140_", "X", row.names(sample_map))
sample_map$off = FALSE
sample_map$off[grepl("N",sample_map$sample.1)] = TRUE
sample_map$off[grepl("M",sample_map$sample.1)] = TRUE
sample_map$off[grepl("Blank",sample_map$sample.1)] = TRUE
sample_map["X1036", "off"] = TRUE


meso_dat = t(data.frame(list(c(1,"shade", "No"),c(2,"shade", "RO"),c(3,"No", "RO"),c(4,"No", "No"),c(5,"shade", "RO"),c(6,"No", "RO"),c(7,"No", "No"),c(8,"shade", "No"),c(9,"No", "No"),c(10,"shade", "RO"),c(11,"No", "RO"),c(12,"No", "No"),c(13,"shade", "RO"),c(14,"shade", "No"),c(15,"No", "No"),c(16,"shade", "No"),c(17,"No", "RO"),c(18,"shade", "RO"),c(19,"shade", "No"),c(20,"No", "RO"),c(21,"ZZ", "ZZ"),c(22,"ZZ", "ZZ"))))

row.names(meso_dat) = NULL
colnames(meso_dat) = c("id", "shade", "RO")


mds = metaMDS(t(data[,row.names(sample_map)[grepl("Exp", sample_map$experiment) & !sample_map$off ]]))$point
mds = cbind(as.data.frame(mds), sample_map[row.names(mds),])


mds$meso = as.numeric(sapply(strsplit(as.character(mds$sample.1), ".", fixed=TRUE), function(x) x[2]))
mds$tpoint = as.numeric(sapply(strsplit(as.character(mds$sample.1), ".", fixed=TRUE), function(x) x[1]))
mds = cbind(mds, meso_dat[mds$meso,])

ggsave("ExpAB.pdf")

mds = metaMDS(t(data[,row.names(sample_map)[ !sample_map$off ]]))$point
mds = cbind(as.data.frame(mds), sample_map[row.names(mds),])

meso2017 =  read.table("../meso2017.csv", sep =",", head = TRUE, row.names = 1)

mds = metaMDS(t(data[,row.names(sample_map)[grepl("2017", sample_map$experiment) & !sample_map$off ]]))$point
mds = cbind(as.data.frame(mds), sample_map[row.names(mds),])
mds = cbind(mds,meso2017[as.numeric(as.character(mds$sample.1)),])

ggplot(mds, aes(MDS1, MDS2, col=as.numeric(nutrients)))+geom_point(size = 5)
