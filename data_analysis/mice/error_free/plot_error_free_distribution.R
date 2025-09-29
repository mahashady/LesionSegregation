library(ggplot2)
library(stringr)
library(gridExtra)
library(ggpubr)
library(ggpmisc)

#file_name <- "../results/fit_normal_mixture_1div_4cells_Egfr_11_14185624_T_A.txt"
#file_name <- "../results/fit_normal_mixture_1div_4cells_Hras_7_145859242_T_A.txt"
#file_name <- "../results/fit_normal_mixture_1div_4cells_Hras_7_145859242_T_C.txt"
file_name <- "../results/fit_normal_mixture_1div_4cells_Braf_6_37548568_A_T.txt"
#file_name <- "../results/fit_beta_mixture_1div_4cells.txt"
#file_name <- "../results/fit_norm_mixture_1div_4cells.txt"
df <- read.table(file_name, header=TRUE, sep="\t")

df <- df[df$n_mut >= 500,]

driver <- str_split_i(str_split_i(file_name, "/", 3), "4cells_", 2)
distribution <- str_split_i(str_split_i(file_name, "/", 3), "_", 2)
print(head(df,10))
print(is.na(df$mu_subclon))
df <- df[!is.na(df$mu_subclon),]
print(head(df))
print(median(df$w_subclon))
setwd("../plots/")
jpeg(filename=paste("error_free_", driver,".", distribution, ".jpeg", sep=""), width=10, height=10, res=300, units='cm')

(ggplot(df, aes(x=w_subclon))
    + geom_histogram(colour = "darkgrey", fill = "white", bins=10)
    + geom_density()
    + geom_vline(xintercept = median(df$w_subclon), col="red", size=0.5)
    + geom_text(x=median(df$w_subclon)+0.1, y= 4, label= paste("median=", round(median(df$w_subclon), 2)),hjust = 0)
    + ggtitle(paste(driver, " (n=", nrow(df),")",sep=""))
    + theme_bw()
    + xlab("proportion of subclonal mutations")
)

dev.off()

setwd("../plots/")
jpeg(filename=paste("mu_subclon_vs_w_subclone_", driver,".", distribution, ".jpeg", sep=""), width=20, height=13, res=300, units='cm')
p1 <- (ggplot(df, aes(x=mu_subclon, y=w_subclon,col=mu_clon/mu_subclon))
    + geom_point()
    + ggtitle(paste(driver, " (n=", nrow(df),")",sep=""))
    + scale_color_gradient2(low="blue", high="red", mid="darkgrey", midpoint=2)
    + theme_bw()
    + theme(legend.position="bottom")
)

p2 <- (ggplot(df, aes(x=mu_clon, y=mu_subclon, col=mu_clon/mu_subclon))
    + geom_point()
    + ggtitle(paste(driver, " (n=", nrow(df),")",sep=""))
    + theme_bw()
#    + ylim(0.2,0.5)
    + stat_poly_line() 
    + scale_color_gradient2(low="blue", high="red", mid="darkgrey", midpoint=2)
    + stat_poly_eq(use_label(c("eq", "R2")))
    + theme(legend.position="bottom")
    + geom_abline(intercept=0, slope=1/2, col="grey")
)
grid.arrange(p2,p1,nrow=1)
dev.off()