library("readxl")
library(reshape2)
library(stringr)
library(viridis)
library(ggplot2)

colors_list=c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
df_hmm_summary = read.table("../results/HMM/Summary.tsv", sep=" ")
colnames(df_hmm_summary) <- c("sample_id", "A_N2T_N_state5_high_vaf", "A_N2T_N_state5_low_vaf", "Median_vaf_autosomes", "Median_vaf_X", "n_sites_as_hmm_state1", "n_sites_as_hmm_state2","n_sites_as_hmm_state3","n_sites_as_hmm_state4","n_sites_as_hmm_state5", "n_sites_multi_hmm_state1", "n_sites_multi_hmm_state2", "n_sites_multi_hmm_state3", "emission_multi_hmm_state1", "emission_multi_hmm_state2", "emission_multi_hmm_state3", "emission_as_hmm_state1", "emission_as_hmm_state2", "emission_as_hmm_state3", "emission_as_hmm_state4", "emission_as_hmm_state5")
df_hmm_summary$sample <- str_split_i(df_hmm_summary$sample_id, "[.]", 1)
head(df_hmm_summary)

df_ploidy_summary = read.table("../results/HMM/Summary_subclone.tsv")
colnames(df_ploidy_summary) <- c("sample", "flexmix_fit", "ploidy")
head(df_ploidy_summary)

c3h_df <- read_excel("../data/41586_2020_2435_MOESM2_ESM.xlsx", sheet = "C3H_tumours")
c3h_df$mice_line = "C3H"
colnames(c3h_df)[1] <- "sample"
head(c3h_df)
cast_df <- read_excel("../data/41586_2020_2435_MOESM2_ESM.xlsx", sheet = "CAST_tumours")
cast_df$mice_line = "CAST"
colnames(cast_df)[1] <- "sample"
head(cast_df)
df_mice_lines_summary = rbind(c3h_df, cast_df)
head(df_mice_lines_summary)

full_table = merge(df_hmm_summary, df_ploidy_summary, by = "sample")
full_table = merge(full_table, df_mice_lines_summary, by = "sample")
nrow(full_table)

full_table$prop_multi_hmm_state1 = full_table$n_sites_multi_hmm_state1/(full_table$n_sites_multi_hmm_state1 + full_table$n_sites_multi_hmm_state2 + full_table$n_sites_multi_hmm_state3)
full_table$prop_multi_hmm_state2 = full_table$n_sites_multi_hmm_state2/(full_table$n_sites_multi_hmm_state1 + full_table$n_sites_multi_hmm_state2 + full_table$n_sites_multi_hmm_state3)
full_table$prop_multi_hmm_state3 = full_table$n_sites_multi_hmm_state3/(full_table$n_sites_multi_hmm_state1 + full_table$n_sites_multi_hmm_state2 + full_table$n_sites_multi_hmm_state3)

full_table$prop_as_hmm_state1 = full_table$n_sites_as_hmm_state1/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state2 = full_table$n_sites_as_hmm_state2/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state3 = full_table$n_sites_as_hmm_state3/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state4 = full_table$n_sites_as_hmm_state4/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state5 = full_table$n_sites_as_hmm_state5/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)

full_table$Pl2<-rep('Dip',length(full_table[,1])) #assign all samples as diploids
full_table$Pl2[(apply(full_table[,c("prop_as_hmm_state1", "prop_as_hmm_state2", "prop_as_hmm_state3", "prop_as_hmm_state4", "prop_as_hmm_state5")] > 0.05, 1, sum)>= 4)]<-'Tetra' #assign tumors that have 4 or more as-HMM states (>0.05 state proportion) as tetraploid

full_table$tetraploids <- 0
full_table[(full_table$ploidy == "tetra" & full_table$flexmix_fit == "Good")|full_table$Pl2 == "Tetra",]$tetraploids <- 1

full_table$Symmetry <- rep('Asymmetric',length(full_table[,1]))
full_table$Symmetry[(apply(full_table[,c("prop_as_hmm_state1", "prop_as_hmm_state5")] < 0.05, 1,sum) >= 2) | full_table$emission_as_hmm_state1 < 0.8 | full_table$emission_as_hmm_state5 > 0.2 | full_table$ploidy == "Symmetric"]<- "Symmetric"

#divisions = full_table[full_table$tetraploids != 1 & full_table$Symmetry != "Symmetric",]
divisions = full_table[full_table$tetraploids != 1 ,]
divisions$division <- "X"
divisions[divisions$prop_multi_hmm_state1 > 0.95,]$division <- "late"
divisions[divisions$prop_multi_hmm_state1 < 0.05,]$division <- "1"
divisions[divisions$Symmetry == "Symmetric",]$division <- "0"
###division

divisions$dist_div2 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(1/4))^2+(as.numeric(x["prop_multi_hmm_state2"])-(2/4))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/4))^2)})
divisions$dist_div3 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(9/16))^2+(as.numeric(x["prop_multi_hmm_state2"])-(6/16))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/16))^2)})
divisions$dist_div4 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(49/64))^2+(as.numeric(x["prop_multi_hmm_state2"])-(14/64))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/64))^2)})
divisions$dist_div5 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(225/256))^2+(as.numeric(x["prop_multi_hmm_state2"])-(30/256))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/256))^2)})
divisions$min_dist_cluser = apply(divisions, 1, function(x) {as.character(which.min(c(as.numeric(x["dist_div2"]),as.numeric(x["dist_div3"]),as.numeric(x["dist_div4"]),as.numeric(x["dist_div5"])))+1)})
divisions$division = ifelse(divisions$division == "X", divisions$min_dist_cluser, divisions$division)
x=divisions[divisions$division != "1" & divisions$division != "0" & !is.na(divisions$division),]
print(head(x[order(-x$prop_multi_hmm_state2),]))

write.csv(divisions, file = "../results/Summary_divisions_with_symmetrical_no_tetra.txt" , row.names = F)

print(table(divisions$division))
#plt_state1_2_clusters_scatter<-(ggplot(divisions[divisions$division != 1 & divisions$division != 0,], aes(x=prop_multi_hmm_state1, y = prop_multi_hmm_state2, col=as.factor(min_dist_cluser))) 
#                                + geom_segment(aes(x = 1/4, y = 0, xend = 1/4, yend = 2/4), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 0, y = 2/4, xend = 1/4, yend = 2/4), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 9/16, y = 0, xend = 9/16, yend = 6/16), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 0, y = 6/16, xend = 9/16, yend = 6/16), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 49/64, y = 0, xend = 49/64, yend = 14/64), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 0, y = 14/64, xend = 49/64, yend = 14/64), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 225/256, y = 0, xend = 225/256, yend = 30/256), linetype="dashed",col="darkgrey")
#                                + geom_segment(aes(x = 0, y = 30/256, xend = 225/256, yend = 30/256), linetype="dashed",col="darkgrey")
#                                + geom_point(x=1/4, y=2/4, col="grey40")
#                                + geom_point(x=9/16, y=6/16, col="grey40")
#                                + geom_point(x=49/64, y=14/64, col="grey40")
#                                + geom_point(x=225/256, y=30/256, col="grey40")
#                                + geom_point(size=1)
#                                + annotate("text", x=0.16, y=0.48, label="II division")
#                                + annotate("text", x=0.46, y=0.35, label="III division")
#                                + annotate("text", x=0.69, y=0.2, label="IV division")
#                                + annotate("text", x=0.81, y=0.1, label="V division")
#                                + theme_bw()
#                                + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#                                + scale_x_continuous(expand = c(0, 0.01), breaks = c(0, 1/4,9/16,49/64,225/256,1.2), labels = c("0", parse(text="(1/2)^2"), parse(text="(3/4)^2"), parse(text="(7/8)^2"), parse(text="(15/16)^2"),"")) 
#                                + scale_y_continuous(expand = c(0, 0.01), breaks = c(0, 2/4,6/16,14/64,30/256,1.2), labels = c("0", parse(text="2*(1/2)*(1/2)"), parse(text="2*(1/4)*(3/4)"), parse(text="2*(1/8)*(7/8)"), parse(text="2*(1/16)*(15/16)"), ""))
#                                + xlab("proportion of sites in state 1")
#                                + ylab("proportion of sites in state 2")
#                                + theme(legend.position = "top")
#                                + scale_color_manual(name="MRCA generation", values=colors_list[3:6])
#)
#ggsave(plt_state1_2_clusters_scatter, file="../plots/multi_hmm_state1_2_scatter_clusters.jpeg", width = 5, height = 5, dpi = 300)


plot_state1_2<-(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.05& divisions$prop_multi_hmm_state2 > 0.05 & divisions$prop_multi_hmm_state3 > 0 ,], aes(x=prop_multi_hmm_state1, y = prop_multi_hmm_state2)) 
                + geom_segment(aes(x = 1/4, y = 0, xend = 1/4, yend = 2/4), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 0, y = 2/4, xend = 1/4, yend = 2/4), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 9/16, y = 0, xend = 9/16, yend = 6/16), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 0, y = 6/16, xend = 9/16, yend = 6/16), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 49/64, y = 0, xend = 49/64, yend = 14/64), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 0, y = 14/64, xend = 49/64, yend = 14/64), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 225/256, y = 0, xend = 225/256, yend = 30/256), linetype="dashed",col="darkgrey")
                + geom_segment(aes(x = 0, y = 30/256, xend = 225/256, yend = 30/256), linetype="dashed",col="darkgrey")
                + geom_point(x=1/4, y=2/4, col="grey40")
                + geom_point(x=9/16, y=6/16, col="grey40")
                + geom_point(x=49/64, y=14/64, col="grey40")
                + geom_point(x=225/256, y=30/256, col="grey40")
                + geom_point(size=1, col="cornflowerblue")
                + annotate("text", x=0.16, y=0.48, label="II division")
                + annotate("text", x=0.46, y=0.35, label="III division")
                + annotate("text", x=0.69, y=0.2, label="IV division")
                + annotate("text", x=0.81, y=0.1, label="V division")
                + theme_bw()
                + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                + scale_x_continuous(expand = c(0, 0.01), breaks = c(0, 1/4,9/16,49/64,225/256,1.2), labels = c("0", parse(text="(1/2)^2"), parse(text="(3/4)^2"), parse(text="(7/8)^2"), parse(text="(15/16)^2"),"")) 
                + scale_y_continuous(expand = c(0, 0.01), breaks = c(0, 2/4,6/16,14/64,30/256,1.2), labels = c("0", parse(text="2*(1/2)*(1/2)"), parse(text="2*(1/4)*(3/4)"), parse(text="2*(1/8)*(7/8)"), parse(text="2*(1/16)*(15/16)"), ""))
)
ggsave(plot_state1_2, file="/home/mandrianova/Lab/Burst_kinetics/mice_pulses_2023/results/Vova_annotation/plots/multi_hmm_state1_2_scatter.jpeg", width = 7, height = 7, dpi = 300)

(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.1,], aes(x=prop_multi_hmm_state1, y = prop_multi_hmm_state3)) 
  + geom_point() 
  + geom_point(x=1/4, y=1/4, col="red")
  + geom_point(x=9/16, y=1/16, col="red")
  + geom_point(x=49/64, y=1/64, col="red")
  + geom_point(x=225/256, y=1/256, col="red")
  + theme_bw())


(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.1,], aes(x=prop_multi_hmm_state1, y = prop_multi_hmm_state3)) 
  + geom_point(col="grey35") 
  + geom_segment(aes(x = 1/4, y = 0, xend = 1/4, yend = 1/4), linetype="dashed",col="blue")
  + geom_segment(aes(x = 0, y = 1/4, xend = 1/4, yend = 1/4), linetype="dashed",col="blue")
  + geom_segment(aes(x = 9/16, y = 0, xend = 9/16, yend = 1/16), linetype="dashed",col="blue")
  + geom_segment(aes(x = 0, y = 1/16, xend = 9/16, yend = 1/16), linetype="dashed",col="blue")
  + geom_point(x=9/16, y=1/16, col="blue")
  + geom_point(x=49/64, y=1/64, col="red")
  + geom_point(x=225/256, y=1/256, col="red")
  + theme_bw())

plot_state1_distr<-(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.1 & divisions$prop_multi_hmm_state2 > 0.05 & divisions$prop_multi_hmm_state3 > 0 & divisions$mice_line=="C3H" ,], aes(x=prop_multi_hmm_state1)) 
                    + geom_histogram(fill="darkgrey")
                    + geom_vline(xintercept = 1/4,linetype="dashed")
                    + geom_vline(xintercept = 9/16, linetype="dashed")
                    + geom_vline(xintercept = 49/64, linetype="dashed")
                    + geom_vline(xintercept = 225/256, linetype="dashed")
                    + xlab("proporion of genome without multi-allelic sites")
                    + ylab("Number of samples")
                    + annotate("text", x=0.30, y=18, label="(1/2)^2",parse=TRUE)
                    + annotate("text", x=0.30, y=16, label="II division")
                    + annotate("text", x=0.61, y=18, label="(3/4)^2",parse=TRUE)
                    + annotate("text", x=0.61, y=16, label="III division")
                    + annotate("text", x=0.81, y=18, label="(7/8)^2",parse=TRUE)
                    + annotate("text", x=0.81, y=16, label="IV division")
                    + annotate("text", x=0.93, y=18, label="(15/16)^2",parse=TRUE)
                    + annotate("text", x=0.93, y=16, label="V division")
                    + theme_bw()
                    + ggtitle("C3H"))
ggsave(plot_state1_distr, file="/home/mandrianova/Lab/Burst_kinetics/mice_pulses_2023/results/Vova_annotation/plots/multi_hmm_state1_distribution_C3H.jpeg", width = 9, height = 5, dpi = 300)


(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.1 & divisions$prop_multi_hmm_state2 > 0.05 & divisions$prop_multi_hmm_state3 > 0,], aes(x=prop_multi_hmm_state2)) 
  + geom_histogram(fill="darkgrey")
  + geom_vline(xintercept = 2/4,linetype="dashed")
  + geom_vline(xintercept = 6/16, linetype="dashed")
  + geom_vline(xintercept = 14/64, linetype="dashed")
  + geom_vline(xintercept = 30/256, linetype="dashed")
  + xlab("proporion of genome without multi-allelic sites")
  + theme_bw()
)

(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.1 & divisions$prop_multi_hmm_state2 > 0.05 & divisions$prop_multi_hmm_state3 > 0,], aes(x=prop_multi_hmm_state3)) 
  + geom_histogram(fill="darkgrey")
  + geom_vline(xintercept = 1/4,linetype="dashed")
  + geom_vline(xintercept = 1/16, linetype="dashed")
  + geom_vline(xintercept = 1/64, linetype="dashed")
  + geom_vline(xintercept = 1/256, linetype="dashed")
  + xlab("proporion of genome without multi-allelic sites")
  + theme_bw()
  + scale_x_continuous(trans="log"))

(ggplot(divisions[divisions$prop_multi_hmm_state1 >= 0.1 & divisions$mice_line == "C3H",], aes(x=prop_multi_hmm_state1)) 
  + geom_histogram(bins=25, fill="grey")
  + geom_vline(xintercept = 2/4,col="red")
  + geom_vline(xintercept = 6/16, col="red")
  + geom_vline(xintercept = 14/64, col="red")
  + geom_vline(xintercept = 30/256, col="red")
  + theme_bw())


result_divisions <- data.frame("Division" = c("0","1","2","3","4","5+"), "C3H" = c(20,108,93,25,29,33), "CAST" = c(0,15,18,4,6,8))
result_divisions_both = data.frame("Division" = c("0","1","2","3","4","5+"), "Both" = result_divisions$C3H + result_divisions$CAST)
result_divisions_both$Both = result_divisions_both$Both/sum(result_divisions_both$Both)
result_divisions$C3H = result_divisions$C3H/sum(result_divisions$C3H)
result_divisions$CAST = result_divisions$CAST/sum(result_divisions$CAST)
result_divisions$Division = as.factor(result_divisions$Division)
result_divisions_melt = melt(result_divisions)
ggplot(result_divisions_melt, aes(x = Division, y=value,fill=variable))+geom_bar(stat = "identity", position = "dodge")
MRCA_distr_by_line <- (ggplot(result_divisions_melt[result_divisions_melt$variable == "CAST"|result_divisions_melt$variable == "C3H",], aes(x = Division, y=value,fill=variable))
                       + geom_bar(stat = "identity", position = "dodge")
                       + theme_bw()
                       + ylab("proportion of samples")
                       + scale_fill_discrete(name="mice line")
                       + xlab("MRCA generation"))
ggsave(MRCA_distr_by_line, file="/home/mandrianova/Lab/Burst_kinetics/mice_pulses_2023/results/Vova_annotation/plots/MRCA_distribution_by_line.jpeg", width = 6, height = 4, dpi = 300)

MRCA_distr_together <- (ggplot(result_divisions_both, aes(x = Division, y=Both, fill=Division))
                        + geom_bar(stat = "identity", position = "dodge", col="#2F4F4F", alpha=0.7)
                        + theme_bw()
                        + ylab("proportion of samples")
                        + xlab("MRCA generation")
                        + scale_fill_manual(values=colors_list)
                        + theme(legend.position = "none"))
ggsave(MRCA_distr_together, file="/home/mandrianova/Lab/Burst_kinetics/mice_pulses_2023/results/Vova_annotation/plots/MRCA_distribution_both_for_poster.jpeg", width = 4, height = 5, dpi = 300)

###########drivers
for_regression <- divisions[!is.na(divisions$knownDrivers),]
for_regression <- for_regression[for_regression$knownDrivers == "Braf"|for_regression$knownDrivers == "Hras"| for_regression$knownDrivers == "Egfr",]
for_regression$division <- ifelse(for_regression$division == "late", 6, for_regression$division)
model = lm(division~ knownDrivers*mice_line, data=for_regression)
