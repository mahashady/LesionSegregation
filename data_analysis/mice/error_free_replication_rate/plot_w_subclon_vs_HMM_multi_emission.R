library(ggplot2)
library(stringr)
library(gridExtra)
library(ggpubr)
library(ggpmisc)

file_subclone <- "../results/fit_norm_mixture_1div_4cells.txt"
df_subclone <- read.table(file_subclone, header=TRUE, sep="\t")
df_subclone <- df_subclone[!is.na(df_subclone$mu_subclon),]
print(nrow(df_subclone))
df_multi <- read.table("../LAD/results/Summary_divisions_with_symmetrical_no_mixtures.txt", header=TRUE, sep=",")
print(nrow(df_multi))

df <- merge(df_subclone, df_multi, by= "sample")
print(nrow(df_multi))


setwd("../plots/")
jpeg(filename="w_subclon_norm_vs_HMMmulti_emission.jpeg", width=30, height=10, res=300, units='cm')

p1 <- (ggplot(df, aes(x=emission_multi_hmm_state3,y=w_subclon, col=mu_clon/mu_subclon))
    + geom_point()
    + stat_poly_line() 
    + stat_poly_eq(use_label(c("eq", "R2")))
    + theme_bw()
    + scale_color_gradient2(low="blue", high="red", mid="darkgrey", midpoint=2)
)

p2 <- (ggplot(df, aes(x=emission_multi_hmm_state3,y=mu_subclon, col=mu_clon/mu_subclon))
    + geom_point()
    + stat_poly_line() 
    + stat_poly_eq(use_label(c("eq", "R2")))
    + theme_bw()
    + scale_color_gradient2(low="blue", high="red", mid="darkgrey", midpoint=2)
)

grid.arrange(p1,p2,nrow=1)
dev.off()
