library(ggplot2)
library(ggpmisc)
library(hash)
library(stringr)

#possible drivers are: 11:14185624_T/A, 6:37548568_A/T, 7:145859242_T/C, 7:145859242_T/A
driver = "7:145859242_T/A" 
compl_nucl <- hash("A"="T", "T"="A" )
opposite_strand <- hash("+"="-", "-"="+")

df_drivers = read.table("../data/Drivers_annotated.csv", header=TRUE, sep=",")
driver_gene = df_drivers[df_drivers$mutID == driver,]$Gene_name
driver_pos = df_drivers[df_drivers$mutID == driver,]$pos
driver_expression = df_drivers[df_drivers$mutID == driver,]$expression_q
driver_strand = df_drivers[df_drivers$mutID == driver,]$Gene_strand
driver_refN = df_drivers[df_drivers$mutID == driver,]$refN
driver_altN = df_drivers[df_drivers$mutID == driver,]$altN

df_MRCA <- read.table("../LAD/results/Summary_divisions_with_symmetrical_no_mixtures.txt", header=TRUE, sep=",")
print(table(df_MRCA$division))

result_samples <- NULL
result_multi <- NULL
result_all <- NULL
result_rates <- NULL
result_division <- NULL
for (sample_name in df_MRCA$sample){
    print(sample_name)
    print(paste("Division=", df_MRCA[df_MRCA$sample==sample_name,]$division, sep=""))
    if (df_MRCA[df_MRCA$sample==sample_name,]$division == 1){
        df_sample <- read.table(paste("../data/mutations_vs_genes_vs_HMM_multi_state/", sample_name, ".gene.HMMmulti_state.nodMat", sep=""), header= TRUE, sep= ",")
        df_sample <- df_sample[df_sample$HMM_multi_state == "A2",]
        multi_sum <- (nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == driver_expression & df_sample$geneStrand == driver_strand & df_sample$refN == driver_refN,])+nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == driver_expression & df_sample$geneStrand == opposite_strand[[driver_strand]] & df_sample$refN == compl_nucl[[driver_refN]],]))
        all_sum <- (nrow(df_sample[df_sample$expression_q == driver_expression & df_sample$geneStrand == driver_strand & df_sample$refN == driver_refN,])+nrow(df_sample[df_sample$expression_q == driver_expression & df_sample$geneStrand == opposite_strand[[driver_strand]] & df_sample$refN == compl_nucl[[driver_refN]],]))
        multi_rate <- multi_sum/all_sum
        if (all_sum > 100){
            result_samples <- c(result_samples, sample_name)
            result_multi <- c(result_multi, multi_sum)
            result_all <- c(result_all, all_sum)
            result_rates <- c(result_rates, multi_rate/2)
            result_division <- c(result_division, df_MRCA[df_MRCA$sample==sample_name,]$division)
        }
    } else if(df_MRCA[df_MRCA$sample==sample_name,]$division == 0){
        next
    } else{
        df_sample <- read.table(paste("../data/mutations_vs_genes_vs_HMM_multi_state/", sample_name, ".gene.HMMmulti_state.nodMat", sep=""), header= TRUE, sep= ",")
        df_sample <- df_sample[df_sample$HMM_multi_state == "A1",]
        multi_sum <- (nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == driver_expression & df_sample$geneStrand == driver_strand & df_sample$refN == driver_refN,])+nrow(df_sample[df_sample$Multi_class == "M" & df_sample$expression_q == driver_expression & df_sample$geneStrand == opposite_strand[[driver_strand]] & df_sample$refN == compl_nucl[[driver_refN]],]))
        all_sum <- (nrow(df_sample[df_sample$expression_q == driver_expression & df_sample$geneStrand == driver_strand & df_sample$refN == driver_refN,])+nrow(df_sample[df_sample$expression_q == driver_expression & df_sample$geneStrand == opposite_strand[[driver_strand]] & df_sample$refN == compl_nucl[[driver_refN]],]))
        multi_rate <- multi_sum/all_sum
        if (all_sum > 100){
            result_samples <- c(result_samples, sample_name)
            result_multi <- c(result_multi, multi_sum)
            result_all <- c(result_all, all_sum)
            result_rates <- c(result_rates, multi_rate)
            result_division <- c(result_division, df_MRCA[df_MRCA$sample==sample_name,]$division)
        }
    }
}

result <- data.frame("sample"=result_samples, "MRCA" = result_division,"n_multi"=result_multi, "n_all"=result_all, "multi_sites_proportion" = result_rates)
result <- result[result$multi_sites_proportion != 0,]
result$division <- ifelse(result$MRCA != "late", result$MRCA, "5")
print(head(result[order(result$division, decreasing=TRUE),]))
setwd("../plots/")
jpeg(filename=paste("repair_rate.", driver_gene, ".", driver_pos, ".", driver_refN, driver_altN, ".jpeg", sep=""), width=15, height=10, res=300, units='cm')

div_median_multi_rate = aggregate(multi_sites_proportion ~ division, data=result, median)
colnames(div_median_multi_rate) <- c("division", "median_multi_rate")
div_median_multi_rate <- div_median_multi_rate[div_median_multi_rate$median_multi_rate != 0,]
div_median_multi_rate$median_multi_rate_log <- log(div_median_multi_rate$median_multi_rate)
m <- lm(log(median_multi_rate) ~ as.numeric(division), data = div_median_multi_rate)
print(summary(m))
slope <- round(coef(m)[[2]],2)
print(slope)
r <- round(1 - exp(slope),2)

(ggplot(result, aes(x=division,y=multi_sites_proportion))
    + geom_jitter(color = "darkgrey", size=1)
    + geom_boxplot(alpha=0.5)
    + theme_bw()
    + ylab("Multiallelic variants rate")
    + scale_y_continuous(trans="log", breaks=c(0.01,0.02, 0.03, 0.04, 0.05, 0.1))
    + stat_poly_line(data=div_median_multi_rate,aes(x=as.numeric(division), y=median_multi_rate),col="black",se = F) 
    + stat_poly_eq(data=div_median_multi_rate,aes(x=as.numeric(division), y=median_multi_rate, label=after_stat(eq.label)))
    + annotate(geom="text", x=4, y = 0.14, label=paste("r=", r, sep=""))
    + ylab("Proportion of multiallelic sites from all variants") 
    + ggtitle(paste(driver_gene, driver, sep="_"))   
)

dev.off()
