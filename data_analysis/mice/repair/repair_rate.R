library(ggplot2)
library(ggpmisc)


df_MRCA <- read.table("../../MRCA/results/Summary_mixed_excluded.txt", header=TRUE, sep=",")
print(table(df_MRCA$division))

result_samples <- NULL
result_rates <- NULL
result_division <- NULL
for (sample_name in df_MRCA$sample){
    print(sample_name)
    print(paste("Division=", df_MRCA[df_MRCA$sample==sample_name,]$division, sep=""))
    if (df_MRCA[df_MRCA$sample==sample_name,]$division == 1){
        df_sample <- read.table(paste("../../data/mutations_vs_genes_vs_HMM_multi_state/", sample_name, ".gene.HMMmulti_state.nodMat", sep=""), header= TRUE, sep= ",")
        df_sample <- df_sample[df_sample$HMM_multi_state == "A2",]
        multi_rate <- nrow(df_sample[df_sample$Multi_class == "M",])/nrow(df_sample)
        result_samples <- c(result_samples, sample_name)
        result_rates <- c(result_rates, multi_rate/2)
        result_division <- c(result_division, df_MRCA[df_MRCA$sample==sample_name,]$division)
    } else if(df_MRCA[df_MRCA$sample==sample_name,]$division == 0){
        next
    } else{
        df_sample <- read.table(paste("../../data/mutations_vs_genes_vs_HMM_multi_state/", sample_name, ".gene.HMMmulti_state.nodMat", sep=""), header= TRUE, sep= ",")
        df_sample <- df_sample[df_sample$HMM_multi_state == "A1",]
        multi_rate <- nrow(df_sample[df_sample$Multi_class == "M",])/nrow(df_sample)
        result_samples <- c(result_samples, sample_name)
        result_rates <- c(result_rates, multi_rate)
        result_division <- c(result_division, df_MRCA[df_MRCA$sample==sample_name,]$division)
    }
}

result <- data.frame("sample"=result_samples, "MRCA" = result_division,"multi_sites_proportion" = result_rates)
result$division <- ifelse(result$MRCA != "late", result$MRCA, "5")
print(head(result))
setwd("../plots/")



jpeg(filename="repair_rate_all.jpeg", width=15, height=10, res=300, units='cm')
div_median_multi_rate = aggregate( multi_sites_proportion ~ division, data=result, median)
m <- lm(log(multi_sites_proportion) ~ as.numeric(division), data = div_median_multi_rate)
print(summary(m))
slope <- coef(m)[[2]]
print(slope)
r <- round(1 - exp(slope),2)
colnames(div_median_multi_rate) <- c("division", "median_multi_rate")
print(head(div_median_multi_rate))
(ggplot(result, aes(x=division,y=multi_sites_proportion))
    + geom_jitter(color = "darkgrey", size=1)
    + geom_boxplot(alpha=0.5)
    + theme_bw()
    + ylab("Multiallelic variants rate")
    + scale_y_continuous(trans="log")
    + stat_poly_line(data=div_median_multi_rate,aes(x=as.numeric(division), y=median_multi_rate),col="black",se = F) 
    + stat_poly_eq(data=div_median_multi_rate,aes(x=as.numeric(division), y=median_multi_rate, label=after_stat(eq.label)))
    + annotate(geom="text", x=4, y = 0.14, label=paste("r=", r, sep=""))
    + ylab("Proportion of multiallelic sites from all variants")    
)

dev.off()
