library(ggplot2)
library(plyr)
library(reshape2)
library(vegan)
library(chron)

data = read.csv("../001_data/ordination for IEG day.csv")
data$dept[data$dept=="Erken"]  = "Limnology program"
data$answer.count = 4-rowSums(apply(data[,4:7], 2, function(x) x==""))
data$why[data$why=="love"]  = "Love"
data$why[data$why=="fun"]  = "Fun"
data$why[data$why=="Don't know"]  = ""

data$obstacle[data$obstacle=="Work-life"]  = "Work-Life"
data$why = factor(data$why)
data$funny = as.factor(data$funny == 1)
data$future[data$future=="post doc"]  = "in Academia"
data$future[data$future=="Reserach"]  = "in Academia"
data$future[data$future=="No idea."]  = ""
data$future[data$future=="Here"]  = "in Academia"
data$future[data$future=="in Academia or far away"]  = "in Academia"
data$date = sapply(strsplit(as.vector(data$Timestamp)," "),'[',1)
data$time = sapply(strsplit(as.vector(data$Timestamp)," "),'[',2)
data$chron = chron(dates = data$date, times = data$time, format=c("d/m/Y", "h:m:s"))
data$job = factor(data$job, levels = levels(data$job)[c(5,2,1,3,4)])

data$dream[data$dream=="Buissnes"]  = "Completly different"
data$dream[data$dream=="A stable reserach position"]  = "An academic career"
data$dream[data$dream=="Im living it"]  = "An academic career"
data$dream[data$dream=="science and socitey"]  = "An academic career"
data$dream[data$dream==""]  = "No idea"
data$dream[data$dream=="teach"]  = "Completly different"
data$dream[data$dream=="to old"]  = "No idea"



means_jobs = ddply(data, .(job), summarize, mean_order=mean(chron, na.rm=TRUE), mean_answers = mean(answer.count, na.rm=TRUE))
ggplot(data, aes(x=chron, fill=job))+ geom_histogram(bins=17)+geom_vline(data=means_jobs, aes(xintercept = mean_order))+facet_wrap(~job, ncol=1)+theme_bw()+ guides(fill=FALSE)+scale_x_chron()+xlab("Registration time")+scale_fill_brewer(palette="Set2")
ggsave("time_by_job.pdf", width=6, height=12)

ggplot(data, aes(x=answer.count, fill=job))+ geom_histogram(bins=5)+geom_vline(data=means_jobs, aes(xintercept = mean_answers))+facet_wrap(~job, ncol=1)+theme_bw()+ guides(fill=FALSE)+scale_fill_brewer(palette="Set1")
ggsave("answers_by_job.pdf", width=6, height=12)

means_dept = ddply(data, .(dept), summarize, mean_order=mean(chron, na.rm=TRUE), mean_answers = mean(answer.count, na.rm=TRUE))
ggplot(data, aes(x=chron, fill=dept))+ geom_histogram(bins=17)+geom_vline(data=means_dept, aes(xintercept = mean_order))+facet_wrap(~dept, ncol=1)+theme_bw()+ guides(fill=FALSE)+scale_x_chron()+xlab("Registration time")+scale_fill_brewer(palette="Set1")
ggsave("time_by_prog.pdf", width=6, height=12)

ggplot(data, aes(x=answer.count, fill=dept))+ geom_histogram(bins=5)+geom_vline(data=means_dept, aes(xintercept = mean_answers))+facet_wrap(~dept, ncol=1)+theme_bw()+ guides(fill=FALSE)
ggsave("answers_by_prog.pdf", width=6, height=12)

good.data = data[,4:7] #[data$answer.count == 4,4:7]

dists = apply(good.data, 1, function(x) apply(good.data, 1, function(y) sum(x != y )/4))
pca = prcomp(dists)

coords = as.data.frame(predict(pca)[,1:4])
coords = cbind(coords, data[row.names(coords),])

#ggplot(coords, aes(x=PC1, y=PC2, col=dream, shape=job))+geom_point(size=5)+theme_bw()
ggplot(coords[as.numeric(coords$job) %in% 2:4 ,], aes(x=PC1, y=PC2, col=dream, shape=job))+geom_point(size=5)+theme_bw()

mds=metaMDS(as.dist(dists))
mds.coords = as.data.frame(mds$points)
mds.coords = cbind(mds.coords, data[row.names(mds.coords),])

ggplot(mds.coords[as.numeric(coords$job) %in% 2:4 ,], aes(x=MDS1, y=MDS2, col=dream, shape=job))+geom_point(size=5)+theme_bw()
ggsave("MDS_dream.pdf", width=10, height=6)
ggplot(mds.coords[as.numeric(coords$job) %in% 2:4 ,], aes(x=MDS1, y=MDS2, col=future, shape=job))+geom_point(size=5)+theme_bw()
ggsave("MDS_future.pdf", width=10, height=6)


whys = melt(sapply(levels(data$job)[2:4], function(x) table(data[data$job == x, "why"])/sum(data$job == x)))
#whys = whys[whys$Var1 != "",]
whys$Var1 = factor(whys$Var1, levels = levels(whys$Var1)[c(2,3,4,5,6,7,1)])
ggplot(whys, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("")+ylab("Proportion")
ggsave("why.pdf", width=10, height=6)

whys = melt(sapply(levels(data$dept), function(x) table(data[data$dept == x, "why"])/sum(data$dept == x)))
#whys = whys[whys$Var1 != "",]

whys$Var1 = factor(whys$Var1)
whys = whys[!is.na(whys$value == whys$value),]

ggplot(whys, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("")+ylab("Proportion")
ggsave("why_dept.pdf", width=10, height=6)



futures = melt(sapply(levels(data$job)[c(2,3,4)], function(x) table(data[data$job == x, "future"])/sum(data$job == x)))
#futures = futures[futures$Var1 != "",]
futures$Var2 = factor(futures$Var2)
futures = futures[futures$value > 0, ]
ggplot(futures, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("")+ylab("Proportion")
ggsave("future_job.pdf", width=10, height=6)

futures = melt(sapply(levels(data$dept)[c(2,4:6)], function(x) table(data[data$dept == x, "future"])/sum(data$dept == x)))
#futures = futures[futures$Var1 != "",]
futures$Var2 = factor(futures$Var2)
futures = futures[!is.na(futures$value) & futures$value > 0, ]
ggplot(futures, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("")+ylab("Proportion")
ggsave("future_dept.pdf", width=10, height=6)



dreams = melt(sapply(levels(data$job), function(x) table(data[data$job == x, "dream"])/sum(data$job == x)))
#dreams = dreams[dreams$Var1 != "",]
#dreams$Var2 = factor(dreams$Var2, levels = levels(dreams$Var2)[c(2,1,3)])
dreams = dreams[dreams$value > 0, ]
ggplot(dreams, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")

dreams = melt(sapply(levels(data$dept), function(x) table(data[data$dept == x, "dream"])/sum(data$dept == x)))
dreams = dreams[dreams$value > 0, ]
ggplot(dreams, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")

funnys = melt(sapply(levels(data$job), function(x) table(data[data$job == x & !is.na(data$funny), "funny"])/sum(data$job[!is.na(data$funny)] == x)))
#funnys = funnys[funnys$Var1 != "",]
funnys$Var2 = factor(funnys$Var2, levels = levels(funnys$Var2))
funnys = funnys[funnys$value > 0, ]
ggplot(funnys, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("Funny")+ylab("Proportion")
ggsave("funny_job.pdf", width=10, height=6)

funnys = melt(sapply(levels(data$dept), function(x) table(data[data$dept == x & !is.na(data$funny), "funny"])/sum(data$dept[!is.na(data$funny)] == x)))
#funnys = funnys[funnys$Var1 != "",]
funnys$Var2 = factor(funnys$Var2)
funnys = funnys[!is.na(funnys$value), ]
ggplot(funnys, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("Funny")+ylab("Proportion")
ggsave("funny_dept.pdf", width=10, height=6)

funnys = melt(sapply(levels(data$why), function(x) table(data[data$why == x & !is.na(data$funny), "funny"])/sum(data$why[!is.na(data$funny)] == x)))
#funnys = funnys[funnys$Var1 != "",]
funnys$Var2 = factor(funnys$Var2)
funnys = funnys[!is.na(funnys$value), ]
ggplot(funnys, aes(x=Var1, fill=Var2, y = value))+geom_histogram(stat = "identity", position="dodge")+theme_bw()+xlab("Funny")+ylab("Proportion")
ggsave("funny_why.pdf", width=10, height=6)
