library(ggplot2)
library(ggpmisc)


df_MRCA <- read.table("../../MRCA/results/Summary_mixed_excluded.txt", header=TRUE, sep=",")
print(table(df_MRCA$division))

result_samples <- NULL
result_rates <- NULL
result_division <- NULL
result_expression <- NULL
for (sample_name in df_MRCA$sample){
    print(sample_name)
    print(paste("Division=", df_MRCA[df_MRCA$sample==sample_name,]$division, sep=""))
    if (df_MRCA[df_MRCA$sample==sample_name,]$division == 1){
        df_sample <- read.table(paste("../../data/mutations_vs_genes_vs_HMM_multi_state/", sample_name, ".gene.HMMmulti_state.nodMat", sep=""), header= TRUE, sep= ",")
        df_sample <- df_sample[df_sample$HMM_multi_state == "A2",]
        multi_rate_q1 <- nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == "q1",])/nrow(df_sample[df_sample$expression_q == "q1",])
        multi_rate_q5 <- nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == "q5",])/nrow(df_sample[df_sample$expression_q == "q5",])
        result_samples <- c(result_samples, sample_name, sample_name)
        result_rates <- c(result_rates, multi_rate_q1/2, multi_rate_q5/2)
        result_division <- c(result_division, df_MRCA[df_MRCA$sample==sample_name,]$division, df_MRCA[df_MRCA$sample==sample_name,]$division)
        result_expression <- c(result_expression, "q1", "q5")
    } else if(df_MRCA[df_MRCA$sample==sample_name,]$division == 0){
        next
    } else{
        df_sample <- read.table(paste("../../data/mutations_vs_genes_vs_HMM_multi_state/", sample_name, ".gene.HMMmulti_state.nodMat", sep=""), header= TRUE, sep= ",")
        df_sample <- df_sample[df_sample$HMM_multi_state == "A1",]
        multi_rate_q1 <- nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == "q1",])/nrow(df_sample[df_sample$expression_q == "q1",])
        multi_rate_q5 <- nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == "q5",])/nrow(df_sample[df_sample$expression_q == "q5",])
        result_samples <- c(result_samples, sample_name, sample_name)
        result_rates <- c(result_rates, multi_rate_q1, multi_rate_q5)
        result_division <- c(result_division, df_MRCA[df_MRCA$sample==sample_name,]$division, df_MRCA[df_MRCA$sample==sample_name,]$division)
        result_expression <- c(result_expression, "q1", "q5")
    }
}

result <- data.frame("sample"=result_samples, "MRCA" = result_division,"multi_sites_proportion" = result_rates,"expression" = result_expression)
result$division <- ifelse(result$MRCA != "late", result$MRCA, "5")
print(head(result))
setwd("../plots/")
jpeg(filename="repair_rate.by_expression.jpeg", width=15, height=10, res=300, units='cm')

#div_median_multi_rate_q1 = aggregate(result[result$expression == "q1",]$multi_sites_proportion, by = list(result[result$expression == "q1",]$division), median)
#colnames(div_median_multi_rate_q1) <- c("division", "median_multi_rate")

#div_median_multi_rate_q5 = aggregate(result[result$expression == "q5",]$multi_sites_proportion, by = list(result[result$expression == "q5",]$division), median)
#colnames(div_median_multi_rate_q5) <- c("division", "median_multi_rate")

div_median_multi_rate = aggregate( multi_sites_proportion ~ division+expression, data=result, median)
print(head(div_median_multi_rate))

m_q1 <- lm(log(multi_sites_proportion) ~ as.numeric(division), data = div_median_multi_rate[div_median_multi_rate$expression == "q1",])
#print(summary(m_q1))
slope_q1 <- coef(m_q1)[[2]]
r_q1 <- 1 - round(exp(slope_q1),2)
print(slope_q1)

m_q5 <- lm(log(multi_sites_proportion) ~ as.numeric(division), data = div_median_multi_rate[div_median_multi_rate$expression == "q5",])
slope_q5 <- coef(m_q5)[[2]]
r_q5 <- round(1 - exp(slope_q5),2)
print(slope_q5)

(ggplot(result, aes(x=division,y=multi_sites_proportion, col=expression))
    + geom_boxplot()
    + geom_point(aes(fill = expression), shape = 21, size = 0.5, position = position_jitterdodge())
    + theme_bw()
    + ylab("Multiallelic variants rate")
    + xlab("division")
    + scale_y_continuous(trans="log", limits = c(0.01, 0.15), breaks=c(0.01,0.02, 0.03, 0.04, 0.05, 0.1))
    + stat_poly_line(data=div_median_multi_rate,aes(x=as.numeric(division), y=multi_sites_proportion),se = F) 
    + stat_poly_eq(data=div_median_multi_rate,aes(x=as.numeric(division), y=multi_sites_proportion, label=after_stat(eq.label)))
    + annotate(geom="text", x=4, y = 0.15, label=paste("r_q1=", r_q1, sep=""))
    + annotate(geom="text", x=4, y = 0.13, label=paste("r_q5=", r_q5, sep=""))
    + ylab("Proportion of multiallelic sites from all variants")    
)

dev.off()
