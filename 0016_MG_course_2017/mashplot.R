library(reshape2)
library(ggplot2)

# load the data
data = read.table("mash_distances.txt", as.is=TRUE)

libraries = levels(as.factor(data$V1))
libraries = sapply(strsplit(libraries,"/"), tail, 1)
libraries = tools::file_path_sans_ext(tools::file_path_sans_ext(libraries))



# extract metadata from library names
#metadata=as.data.frame(t(as.data.frame(strsplit(sub(".gz","",sub(".fastq","",sub(".fq","",libraries))),"_"), optional=FALSE)))

#make row.names and colnames
#row.names(metadata) = libraries
#colnames(metadata) = c("site","day","read")
#metadata$day = as.numeric(metadata$day)

# turn the data into a matrix
pdata = data[,1:3]
distance_matrix = dcast(pdata, V1~V2)
row.names(distance_matrix) = distance_matrix$V1
distance_matrix = distance_matrix[,-1]

# compute a PCA
pca = data.frame(predict(prcomp(distance_matrix))[,1:2])
pca$Row.names = row.names(pca)
#pca = merge(pca,metadata, by.x = 'row.names', by.y = 'row.names')

# make pretty plots
#p1 = ggplot(pca, aes(x=PC1, y=PC2, shape=site, col=day, label=Row.names))+scale_color_gradient(low="purple",high="green") + geom_label()
p1 = ggplot(pca, aes(x=PC1, y=PC2, label=Row.names)) + geom_label()

ggsave("PCA_of_Mash_distances.pdf", p1)
