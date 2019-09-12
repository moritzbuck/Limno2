library(ggplot2)
library(chron)

ysi = read.csv("YSI.csv", as.is=TRUE)
ysi$chron = chron(dates = ysi$Date, times = ysi$Time, format=c("m/d/y", "h:m:s"))
ysi = ysi[ysi$meso_ID != "test",]
treats = read.csv("meso_treats.csv")

ysi = merge(ysi, treats, by="meso_ID")
ggplot(ysi, aes(x=as.POSIXct(chron), y=Temp..C, col=meso_ID, group=depth))+geom_line()+geom_point()+facet_grid(treatment~nutrients)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")
ggsave("Temperature.pdf", width = 12, height = 8)
ggplot(ysi, aes(x=as.POSIXct(chron), y=ODO.mg.L, col=as.factor(depth), group=depth))+geom_line()+geom_point()+facet_grid(treatment~nutrients)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")
ggsave("Dissolved_Oxygen.pdf", width = 12, height = 8)
ggplot(ysi, aes(x=as.POSIXct(chron), y=Chlorophyll.µg.L, col=as.factor(depth), group=depth))+geom_line()+geom_point()+facet_grid(treatment~nutrients)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_y_log10()+xlab("")
ggsave("Dissolved_chlorophyll.pdf", width = 12, height = 8)
ggplot(ysi, aes(x=as.POSIXct(chron), y=ORP.mV, col=as.factor(depth), group=depth))+geom_line()+geom_point()+facet_grid(treatment~nutrients)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")
ggsave("ReductivePotential.pdf", width = 12, height = 8)
ggplot(ysi, aes(x=as.POSIXct(chron), y=Turbidity.FNU, col=as.factor(depth), group=depth))+geom_line()+geom_point()+facet_grid(treatment~nutrients)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")+scale_y_log10()
ggsave("Turbidity.pdf", width = 12, height = 8)
