library(ggplot2)
library(gridExtra)

df_beta <- read.table("/workspace/projects/lesion_segregation/mice/error_free/results/fit_beta_mixture_1div_4cells.txt", header=TRUE, sep="\t")
df_norm <- read.table("/workspace/projects/lesion_segregation/mice/error_free/results/fit_norm_mixture_1div_4cells.txt", header=TRUE, sep="\t")

df <- merge(df_beta, df_norm, by="sample", suffixes=c(".beta", ".norm"))
print(head(df))

p1 <- (ggplot(df,aes(x=mu_subclon.beta, y=mu_subclon.norm))
    + geom_point()
    + theme_bw()
)

p2 <- (ggplot(df,aes(x=mu_clon.beta, y=mu_clon.norm))
    + geom_point()
    + theme_bw()
)

p3 <- (ggplot(df, aes(x=w_subclon.beta, y=w_subclon.norm))
    + geom_point()
    + theme_bw()
)

p4 <- (ggplot(df, aes(x=w_clon.beta, y=w_clon.norm))
    + geom_point()
    + theme_bw()
)
setwd("../plots/")
jpeg(filename="beta_vs_norm_fitting.jpeg", width=30, height=15, res=300, units='cm')
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()
