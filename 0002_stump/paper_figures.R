ggplot(mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")+ggtitle("All 16s")

ggplot(pot_mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")+ggtitle("Potential methylator 16s")


ggplot(sub_mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")+ggtitle("HgcA samples 16s")

ggplot(hgca_mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")+ggtitle("All HgcA")

ggplot(methyl_frac, aes(x=Name, col=Name, y =fraction))+geom_boxplot()+scale_y_log10()+geom_jitter(width = 0.1)+ggtitle("proportion of potential methylators")


clean_style = theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank())

ggplot(family_shaped, aes(x=Var2,y=value, fill=Var1))+geom_col()+scale_fill_manual(values=cols)+guides(fill=guide_legend(ncol=1))+facet_grid(~Name, scale="free_x")+guides(fill=guide_legend(ncol=3))+ ggtitle("Families by samples")+clean_style
ggsave("families.pdf", width = 14, height=8)

ggplot(family_shaped[family_shaped$pot_methy, ], aes(x=Var2,y=value, fill=Var1))+geom_col()+scale_fill_manual(values=cols)+guides(fill=guide_legend(ncol=1))+facet_grid(~Name, scale="free_x")

ggplot(hgca_shaped, aes(x=variable,y=value, fill=taxa))+geom_col()+facet_grid(~Name,scales="free_x")+theme_bw()+scale_fill_manual(values=cols2)
