ggplot(mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")

ggplot(pot_mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")


ggplot(sub_mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")

ggplot(hgca_mds_points, aes(x=MDS1, y=MDS2, col=Ratio.MeHg, shape=Name))+geom_point(size=3)+theme_bw()+scale_color_continuous(trans="log10")


ggplot(methyl_frac, aes(x=Name, col=Name, y =fraction))+geom_boxplot()+scale_y_log10()+geom_jitter(width = 0.1)
