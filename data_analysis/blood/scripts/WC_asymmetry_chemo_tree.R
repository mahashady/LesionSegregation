library(ape)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
source("/home/mandrianova/Lab/Burst_kinetics/Mike/scripts/plot_tree.R")
source("/home/mandrianova/Lab/Burst_kinetics/Mike/scripts/tree_functions.R")
#data_dir = "/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/from_Mendeley_last_version/input_data/input_data/EM/"
data_dir = "/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/from_Mendeley_last_version/input_data/input_data/EM/"
ID="PX001_2_01"
#Load tree
tree_file_path = paste0(data_dir,"tree_",ID,"_standard_rho01.tree")
tree <- di2multi(read.tree(tree_file_path))
load("/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/from_Mendeley_last_version/input_data/input_data/EM/annotated_mut_set_PX001_2_01_standard_rho01")
mutations <- filtered_muts$COMB_mats.tree.build$mat
NV <- filtered_muts$COMB_mats.tree.build$NV
#NR <- filtered_muts$COMB_mats.tree.build$NR

clade_WC_as_by_chrom_based_on_annotation<- function(clade_node, mutations) {
  mutations_shared <- mutations[mutations$node == clade_node,]
  mutations_shared$chrom <- str_split_i(mutations_shared$mut_ref, "-", 1)
  mutations_shared$ref <- str_split_i(mutations_shared$mut_ref, "-", 3)
  mutations_shared$alt <- str_split_i(mutations_shared$mut_ref, "-", 4)
  mutations_shared$substitution <- paste(str_split_i(mutations_shared$mut_ref, "-", 3), str_split_i(mutations_shared$mut_ref, "-", 4), sep=">")
  #subclade_1_result <- as.data.frame(mutations_shared %>% group_by(chrom, substitution) %>%   summarise(number = n()))
  subclade_1_result <- as.data.frame(mutations_shared %>% group_by(chrom, Ref) %>%   summarise(number = n()))
  print(sum(subclade_1_result$number))
  #subclade_1_result <- subclade_1_result[subclade_1_result$substitution == "T>A" |subclade_1_result$substitution == "A>T",]
  subclade_1_result <- subclade_1_result[subclade_1_result$Ref == "T" |subclade_1_result$Ref == "A",]
  subclade_1_result$node = clade_node
  return(subclade_1_result)
}

clade_WC_as_by_chrom <- function(clade_node, NV) {
  clade_colonies <- tree$tip.label[get_all_node_leaves(clade_node,tree)] 
  NV <- as.data.frame(NV)
  NV_pos <- as.data.frame(NV[,c(clade_colonies)])
  NV_neg <- as.data.frame(NV[,-which(names(NV) %in% clade_colonies)])
  NV_1_shared <- NV[apply(NV_pos, 1, function(row) all(row >= 1)) & apply(NV_neg, 1, function(row) all(row <=2)), ]
  NV_1_shared$chrom <- str_split_i(rownames(NV_1_shared), "-", 1)
  NV_1_shared$ref <- str_split_i(rownames(NV_1_shared), "-", 3)
  NV_1_shared$alt <- str_split_i(rownames(NV_1_shared), "-", 4)
  NV_1_shared$substitution <- paste(str_split_i(rownames(NV_1_shared), "-", 3), str_split_i(rownames(NV_1_shared), "-", 4), sep=">")
  subclade_1_result <- as.data.frame(NV_1_shared %>% group_by(chrom, substitution) %>%   summarise(number = n()))
  print(sum(subclade_1_result$number))
  subclade_1_result <- subclade_1_result[subclade_1_result$substitution == "T>A" |subclade_1_result$substitution == "A>T",]
  subclade_1_result$node = clade_node
  return(subclade_1_result)
}


setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
jpeg(filename="WC_as_clades_T>N_A>N.jpg", width=17, height=15, res=300, units='cm')
(ggplot(as_result_for_plot, aes(x=number.T, y=number.A, col=as.factor(node)))
  + geom_point()
  + facet_wrap(~chrom)
  + geom_abline(intercept=0, slope=1)
  + theme_bw()

)

dev.off()


#plot WCas oer node
as_result <- NULL
clades_of_interest <- c(418, 420, 419,432, 447)
for (clade in clades_of_interest){
  print(clade)
  as_clade <- clade_WC_as_by_chrom_based_on_annotation(clade, mutations)
  as_clade$node=clade
  as_result <- rbind(as_result, as_clade)
}
as_result_T <- as_result[as_result$Ref == "T",]
as_result_A <- as_result[as_result$Ref == "A",]
as_result_for_plot <- merge(as_result_T, as_result_A, by=c("chrom", "node"), suffixes = c(".T", ".A"))
as_result_for_plot$node <- as.factor(as_result_for_plot$node)  


setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
jpeg(filename="WC_as_internal_clades_based on_annotation_T>N_A>N.jpg", width=17, height=15, res=300, units='cm')

(ggplot(as_result_for_plot, aes(x=number.T, y=number.A, col=chrom))
  + geom_point()
  + facet_wrap(~node)
  + geom_abline(intercept=0, slope=1)
  + theme_bw()
)
dev.off()

as_result <- NULL
clades_of_interest <- c(420, 419,432, 447)
for (clade in clades_of_interest){
  print(clade)
  as_clade <- clade_WC_as_by_chrom(clade, NV)
  as_clade$node=clade
  as_result <- rbind(as_result, as_clade)
}
as_result_TA <- as_result[as_result$substitution == "T>A",]
as_result_AT <- as_result[as_result$substitution == "A>T",]
as_result_for_plot <- merge(as_result_TA, as_result_AT, by=c("chrom", "node"), suffixes = c(".TA", ".AT"))
as_result_for_plot$node <- as.factor(as_result_for_plot$node)  


setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
jpeg(filename="WC_as_internal_clades.jpg", width=17, height=15, res=300, units='cm')
(ggplot(as_result_for_plot, aes(x=number.T, y=number.A))
  + geom_point()
  + facet_wrap(~node)
  + geom_abline(intercept=0, slope=1)
  + theme_bw()
)
dev.off()

as_result_for_plot$WC_as <- as_result_for_plot$number.T/as_result_for_plot$number.A


(ggplot(as_result_for_plot, aes(x=as.factor(Chromosome), y=WC_as, fill=node))
  + geom_bar(stat="identity", position="dodge")
  #+ geom_point()
  + scale_y_continuous(trans="log", breaks=c(0.5,2))
  + theme_bw()
  + geom_hline(yintercept=1)
)



#### all pair asymmetries

clades_of_interest <- c(420,432)
#clades_of_interest <- c(418,420, 419,432, 447)
#clades_of_interest <- c(267,268, 304)
#clades_of_interest <- c(319,318, 327)
#clades_of_interest <- c(378,377, 380)
#clades_of_interest <- c(268,269, 291)
all_pairs <- as.data.frame(t(combn(clades_of_interest, 2)))
#as_first_clade <- clade_WC_as_by_chrom(first_clade_node, NV)
#as_second_clade <- clade_WC_as_by_chrom(second_clade_node, NV)
p <- list()

for (i in 1:nrow(all_pairs)){
  first_clade_node = all_pairs[i,1]
  second_clade_node = all_pairs[i,2]
  print(paste0(first_clade_node,"-", second_clade_node))
  as_first_clade <- clade_WC_as_by_chrom_based_on_annotation(first_clade_node, mutations)
  as_second_clade <- clade_WC_as_by_chrom_based_on_annotation(second_clade_node, mutations)
  as_result <- rbind(as_first_clade, as_second_clade)
  
  #as_result_TA <- as_result[as_result$substitution == "T>A",]
  #as_result_AT <- as_result[as_result$substitution == "A>T",]
  #as_result_for_plot <- merge(as_result_TA, as_result_AT, by=c("chrom", "node"), suffixes = c(".TA", ".AT"))
  #as_result_for_plot$WC_as <- as_result_for_plot$number.TA/as_result_for_plot$number.AT
  as_result_T <- as_result[as_result$Ref == "T",]
  as_result_A <- as_result[as_result$Ref == "A",]
  as_result_for_plot <- merge(as_result_T, as_result_A, by=c("chrom", "node"), suffixes = c(".T", ".A"))
  
  as_result_for_plot$WC_as <- as_result_for_plot$number.T/as_result_for_plot$number.A
  
  #as_result_for_plot <- as_result_for_plot[as_result_for_plot$number.TA >5 & as_result_for_plot$substitution.AT > 5,]
  #as_result_for_plot$node <- as.factor(as_result_for_plot$node)
  
  
  as_result_node1 = as_result_for_plot[as_result_for_plot$node == first_clade_node,]
  as_result_node2 = as_result_for_plot[as_result_for_plot$node == second_clade_node,]
  as_result_for_plot_between_nodes <- merge(as_result_node1,as_result_node2, by="chrom", suffixes = c(".x", ".y"))
  as_result_for_plot_between_nodes$Total_mut = rowSums(as_result_for_plot_between_nodes[,grepl("number",colnames(as_result_for_plot_between_nodes))])
  #as_result_for_plot_between_nodes$used_for_pval <- ifelse(as_result_for_plot_between_nodes$Total_mut < 50, "no", "yes")
  as_result_for_plot_between_nodes$used_for_pval <- "yes"
  #as_result_for_plot_between_nodes[(as_result_for_plot_between_nodes$WC_as.x < 1.5 & as_result_for_plot_between_nodes$WC_as.x > 2/3) | (as_result_for_plot_between_nodes$WC_as.y < 1.5 & as_result_for_plot_between_nodes$WC_as.y > 2/3) ,]$used_for_pval <- "no"
  print((as_result_for_plot_between_nodes))
  #df_for_pval <- as_result_for_plot_between_nodes[((as_result_for_plot_between_nodes$WC_as.x > 1.5 & as_result_for_plot_between_nodes$WC_as.y > 1.5) | (as_result_for_plot_between_nodes$WC_as.x > 1.5 & as_result_for_plot_between_nodes$WC_as.y < 2/3) | (as_result_for_plot_between_nodes$WC_as.x < 2/3 & as_result_for_plot_between_nodes$WC_as.y > 1.5) | (as_result_for_plot_between_nodes$WC_as.x < 2/3 & as_result_for_plot_between_nodes$WC_as.y < 2/3)),]
  print(nrow(df_for_pval))
  
  df_for_pval$WC_as.x.bin <- ifelse(df_for_pval$WC_as.x > 1, "T", "A")
  df_for_pval$WC_as.y.bin <- ifelse(df_for_pval$WC_as.y > 1, "T", "A")
  df_for_pval$WC_as.x.y.bin <- paste0(df_for_pval$WC_as.x.bin,df_for_pval$WC_as.y.bin )
  print(table(df_for_pval$WC_as.x.y.bin))
  counts_quadro <- table(df_for_pval$WC_as.x.y.bin)
  n_AA = ifelse(!is.na(counts_quadro["AA"]), counts_quadro["AA"], 0)
  n_AT = ifelse(!is.na(counts_quadro["AT"]), counts_quadro["AT"], 0)
  n_TA = ifelse(!is.na(counts_quadro["TA"]), counts_quadro["TA"], 0)
  n_TT = ifelse(!is.na(counts_quadro["TT"]), counts_quadro["TT"], 0)
  df_for_pval_agg <- matrix(c(n_AA, n_AT, n_TA, n_TT), nrow=2)
  print(df_for_pval_agg)
  pval_fisher_cor <- round(fisher.test(df_for_pval_agg,alternative = "greater")$p.value,3)
  pval_fisher_anti <- round(fisher.test(df_for_pval_agg,alternative = "less")$p.value,3)
  print(pval_fisher_cor)
  print(pval_fisher_anti)
  pval_sperman <- round(cor.test(as_result_for_plot_between_nodes$WC_as.x, as_result_for_plot_between_nodes$WC_as.y, method=c("spearman"))$p.value,3)
  as_result_for_plot_between_nodes$used_for_pval <- factor(as_result_for_plot_between_nodes$used_for_pval,levels=c("yes", "no"))
  color_map <- c("yes" = "black", "no" = "lightgrey")
  plot_pair <- (ggplot(as_result_for_plot_between_nodes, aes(x=WC_as.x, y=WC_as.y, size=Total_mut))
    + geom_point()
    + geom_hline(yintercept=1)
    + geom_vline(xintercept=1)
    + theme_bw()
    + scale_x_continuous(trans="log", breaks=c(0.5,2), limits=c(1/15, 15))
    + scale_y_continuous(trans="log", breaks=c(0.5,2), limits=c(1/15, 15))
    + xlab(paste0("WC node", first_clade_node))
    + ylab(paste0("WC node", second_clade_node))
    #+  annotate("text", x=1/15, y=25/2, label= paste0("pval_cor=",pval_fisher_cor), hjust=0)
    #+  annotate("text", x=8/15, y=20/2, label= paste0("pval_fisher=",pval_fisher_anti), hjust=0, size=3)
    +  annotate("text", x=2, y=15/2, label= paste0("pval=",pval_sperman), hjust=0, size=3)
    + scale_color_manual(values=color_map)
  )
  p[[i]] <- plot_pair
  
  
}
do.call(grid.arrange,p)

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
png("WC_plots_420_432.png", width = 10, height = 8, units = "cm", res = 300)
#png("WC_plots_378_377_380.png", width = 30, height = 8, units = "cm", res = 300)
do.call(grid.arrange, c(p, nrow = 1))
dev.off()
