library(phytools)
library(phyloseq)
library(stringr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggrepel)
library(rstatix)
library(ggpubr)


source("/home/mandrianova/Lab/Burst_kinetics/Mike/scripts/tree_functions.R")
source("/home/mandrianova/Lab/Burst_kinetics/Mike/scripts/simulate_tree_controls.R")

data_dir = "/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/trees_and_mutations/"
df_muts <- read.table("/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/trees_and_mutations/mutations_filtered.tsv", sep = "\t", header = T)
#hspc_samples_table = c("PD4781","PD5117","PD5163","PD5179","PD5182","PD5847","PD6629", "PD6646","PD7271","PD9478","PD34493_FINAL","PD41276_v2_FINAL_rho05","PD41305_v3_FINAL","AX001_4_01","CB001_3_01","CB002_2_01","KX001_4_01","KX002_2_01","KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01","KX009_1_01","KX010_1_01","SX001_5_01")
#hspc_samples = c("PD4781","PD5117","PD5163","PD5179","PD5182","PD5847","PD6629", "PD6646","PD7271","PD9478","PD34493_FINAL_noMixed","PD41276_v2_FINAL_rho05_noMixed","PD41305_v3_FINAL_noMixed","AX001_4_01_standard_rho01","CB001_3_01_standard_rho01","CB002_2_01_standard_rho01","KX001_4_01_standard_rho01","KX002_2_01_standard_rho01","KX003_5_01_standard_rho01","KX004_5_01_standard_rho01","KX007_2_01_standard_rho01","KX008_2_01_standard_rho01","KX009_1_01_standard_rho01","KX010_1_01_standard_rho01","SX001_5_01_standard_rho01")
hspc_samples_table = c("PD34493_FINAL","PD41276_v2_FINAL_rho05","PD41305_v3_FINAL","AX001_4_01","CB001_3_01","CB002_2_01","KX001_4_01","KX002_2_01","KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01","KX009_1_01","KX010_1_01","SX001_5_01")
hspc_samples = c("PD34493_FINAL_noMixed","PD41276_v2_FINAL_rho05_noMixed","PD41305_v3_FINAL_noMixed","AX001_4_01_standard_rho01","CB001_3_01_standard_rho01","CB002_2_01_standard_rho01","KX001_4_01_standard_rho01","KX002_2_01_standard_rho01","KX003_5_01_standard_rho01","KX004_5_01_standard_rho01","KX007_2_01_standard_rho01","KX008_2_01_standard_rho01","KX009_1_01_standard_rho01","KX010_1_01_standard_rho01","SX001_5_01_standard_rho01")
hspc_ids = c("PD34493","PD41276","PD41305","AX001_4_01","CB001_3_01","CB002_2_01","KX001_4_01","KX002_2_01","KX003_5_01","KX004_5_01","KX007_2_01","KX008_2_01","KX009_1_01","KX010_1_01","SX001_5_01")
drivers = read.table("/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/driver_table.csv", header = T, sep=",")


########## Test for lesion survival distances
all_dist <- NULL
lesion_durations <- NULL
all_controls_duraion<-NULL
for (sample in hspc_samples){
  print(sample)
  tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  table_sample = str_split_i(sample, "_standard",1)
  table_sample = str_split_i(table_sample, "_noMixed",1)
  print(table_sample)
  sample_pvvs = df_muts[df_muts$Sample_ID == table_sample,]
  sample_pvvs <- sample_pvvs[sample_pvvs$Type =="PVV" & sample_pvvs$Class == "PASS",]
  print(head(sample_pvvs))
  lesion_durations <- c(lesion_durations, sample_pvvs$lesion_duration)
  print(nrow(sample_pvvs))
  if (nrow(sample_pvvs) != 0){
    end_nodes = sample_pvvs$lesion_repair_node
    print(end_nodes)
    start_nodes = sample_pvvs$lesion_node
    for (i in 1:length(end_nodes)){
      interdist = 0 # the number of nodes between lesion_node and lesion_repair_node
      print(start_nodes[i])
      print(end_nodes[i])
      cur_node = end_nodes[i]
      parent_node = 0
      while (parent_node!= start_nodes[i]){
        parent_node = getAncestors(tree, cur_node, type="parent")
        interdist <- interdist + 1 
        cur_node <- parent_node
      }
      print("Lesion finished")
      print(interdist)
      if (nodeheight(tree, cur_node) > 60){
        all_dist <- c(all_dist, interdist)
        print("_______________")
        print("Observed interdist")
        print(interdist)
        print("Observed duration")
        print(sample_pvvs$lesion_duration[i])
        controls <- simulate_controls_duration(tree, end_nodes[i], interdist, 10)
        print("Simulated durations")
        print(controls)
        all_controls_duraion <- c(all_controls_duraion, controls)
      }
      
    }
  }
}
print(all_dist)
df_lesion_duration = data.frame("value"=lesion_durations, type = "lesion")
df_control_duration = data.frame("value"=all_controls_duraion, type = "control")
df_result_duration = rbind(df_lesion_duration, df_control_duration)
control_mean_duration=mean(df_result_duration[df_result_duration$type=="control",]$value)
lesion_mean_duration=mean(df_result_duration[df_result_duration$type=="lesion",]$value)

p_violin_duration<-(ggplot(df_result_duration, aes(y=value, x=type, fill=type)) 
                  + geom_violin() 
                  + geom_jitter(size=0.5,col="gray35")
                  + theme_bw()
                  + ylab("Time in number of mutations")
                  + scale_y_continuous(trans="log",breaks = c(1,10,100))
                  + scale_fill_manual(values=c("grey","coral3"))
                  + theme(legend.position = "none"))


p_hist_duration<-(ggplot(df_result_duration, aes(x=value,fill=type)) 
                + geom_density(alpha=0.8) 
                + theme_bw()
                + xlab("Time in number of mutations")
                + scale_fill_manual(values=c("darkgrey","coral3"), name="")
                + geom_vline(xintercept = lesion_mean_duration, col="coral3", linetype="dashed")
                + geom_vline(xintercept = control_mean_duration, col="darkgrey", linetype="dashed")
)
ks_duration = ks.test(df_result_duration[df_result_duration$type=="lesion",]$value,df_result_duration[df_result_duration$type=="control",]$value)
pval_duration = format(ks_duration$p.value,digits=3,scientific = T)

p_cumul_duration<- (ggplot(df_result_duration, aes(x=value, color=type, y = 1 - after_stat(y))) 
                  + stat_ecdf(geom="step") 
                  + theme_bw()
                  + xlab("Time in number of mutations")
                  + scale_color_manual(values=c("darkgrey","coral3"))
                  + ylab("1-ecdf")
                  + annotate(geom="text",x=500, y=0.5, label=paste("pval=",pval_duration))
                  + theme(legend.position = "none")
)

########## Test for number of node leafs #######################
all_children_numbers <- NULL
lesion_nodes <- NULL
lesion_repair_nodes <- NULL
all_controls<-NULL
for (sample in hspc_samples){
  print(sample)
  tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  table_sample = str_split_i(sample, "_standard",1)
  table_sample = str_split_i(table_sample, "_noMixed",1)
  print(table_sample)
  sample_pvvs = df_muts[df_muts$Sample_ID == table_sample,]
  sample_pvvs <- sample_pvvs[sample_pvvs$Type =="PVV" & sample_pvvs$Class == "PASS",]
  print(head(sample_pvvs))
  lesion_nodes <- c(lesion_nodes, sample_pvvs$lesion_node)
  lesion_repair_nodes <- c(lesion_repair_nodes, sample_pvvs$lesion_repair_node)
  print(nrow(sample_pvvs))
  if (nrow(sample_pvvs) != 0){
    end_nodes = sample_pvvs$lesion_repair_node
    print(end_nodes)
    start_nodes = sample_pvvs$lesion_node
    for (i in 1:length(end_nodes)){
      cur_node = end_nodes[i]
      N_leaves = length(get_all_node_leaves(cur_node, tree))
      print("Observed")
      print(N_leaves)
      controls <- simulate_controls_leaves(tree, cur_node, 50, 5)
      #add observed value only if we were able to find the corresponding control
      if (length(controls) != 0 & nodeheight(tree,cur_node)>60){      
        print("Simulated leaves")
        print(controls)
        all_children_numbers <- c(all_children_numbers, N_leaves)
        all_controls <- c(all_controls, controls)
      }
    }
  }
}
df_lesion_leaves = data.frame("value"=all_children_numbers, type = "lesion")
df_control_leaves = data.frame("value"=all_controls, type = "control")
print(nrow(df_lesion_leaves))
print(nrow(df_control_leaves))
df_result_leaves = rbind(df_control_leaves, df_lesion_leaves)
control_mean_leaves = mean(df_result_leaves[df_result_leaves$type=="control",]$value)
lesion_mean_leaves = mean(df_result_leaves[df_result_leaves$type=="lesion",]$value)

p_violin_leaves<-(ggplot(df_result_leaves, aes(y=value, x=type, fill=type)) 
  + geom_violin() 
  + geom_jitter(size=0.5,col="gray35")
  + theme_bw()
  + ylab("Number of leaves")
  + scale_y_continuous(trans="log",breaks = c(1,10,100))
  + scale_fill_manual(values=c("grey","coral3"))
  + theme(legend.position = "none"))

result_counts_leaves = as.data.frame(table(df_result_leaves))
result_counts_leaves$Frequency = ifelse(result_counts_leaves$type=="control", result_counts_leaves$Freq/sum(result_counts_leaves[result_counts_leaves$type=="control",]$Freq), result_counts_leaves$Freq/sum(result_counts_leaves[result_counts_leaves$type=="lesion",]$Freq))
p_hist_leaves<-(ggplot(result_counts_leaves, aes(x=value, y=Frequency, fill=type)) 
  + geom_histogram(position="dodge", stat="identity") 
  + theme_bw()
  + xlab("Number of leaves")
  + scale_fill_manual(values=c("darkgrey","coral3"))
  + geom_vline(xintercept = lesion_mean_leaves, col="coral3", linetype="dashed")
  + geom_vline(xintercept = control_mean_leaves, col="darkgrey", linetype="dashed")
)

ks_leaves = ks.test(df_result_leaves[df_result_leaves$type=="lesion",]$value,df_result_leaves[df_result_leaves$type=="control",]$value)
pval_leaves = format(ks_leaves$p.value,digits=3,scientific = T)

p_cumul_leaves<- (ggplot(df_result_leaves, aes(x=value, color=type, y = 1 - after_stat(y))) 
  + stat_ecdf(geom="step") 
  + theme_bw()
  + xlab("Number of leaves")
  + scale_color_manual(values=c("darkgrey","coral3"))
  + ylab("1-ecdf")
  + annotate(geom="text",x=50, y=0.5, label=paste("pval=",pval_leaves))
)
########## Test for number of branching points #######################
time_after = 200
#regime could be start or end depending if we analyze origin or repair node of lesion
regime = "end"
number_branching <- NULL
lesion_nodes <- NULL
lesion_repair_nodes <- NULL
all_controls_branching<-NULL
for (sample in hspc_samples){
  print(sample)
  tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  table_sample = str_split_i(sample, "_standard",1)
  table_sample = str_split_i(table_sample, "_noMixed",1)
  print(table_sample)
  sample_pvvs = df_muts[df_muts$Sample_ID == table_sample,]
  sample_pvvs <- sample_pvvs[sample_pvvs$Type =="PVV" & sample_pvvs$Class == "PASS",]
  print(head(sample_pvvs))
  lesion_nodes <- c(lesion_nodes, sample_pvvs$lesion_node)
  lesion_repair_nodes <- c(lesion_repair_nodes, sample_pvvs$lesion_repair_node)
  print(nrow(sample_pvvs))
  if (nrow(sample_pvvs) != 0){
    end_nodes = sample_pvvs$lesion_repair_node
    print(end_nodes)
    start_nodes = sample_pvvs$lesion_node
    for (i in 1:length(end_nodes)){
      if (regime == "start"){
        cur_node = start_nodes[i]
      }else{
        cur_node = end_nodes[i]
        }
      check_longest_offspring = FALSE
      all_cur_node_leaves = get_all_node_leaves(cur_node, tree)
      for (leaf in all_cur_node_leaves){
        if (nodeheight(tree,leaf)-nodeheight(tree,cur_node)>time_after){
          check_longest_offspring = TRUE #this is to check that we have enough time after. If not we will have less branchngs by definition
        }
      }
      if (check_longest_offspring & nodeheight(tree, cur_node)>60){
        N_branchings = length(get_branching_points_in_time(tree, cur_node, time_after))
        print("Observed")
        print(N_branchings)
        if (regime == "start"){
          controls_branching <- simulate_controls_branching_start_node(tree, cur_node, time_after, 5)
        }else{
          controls_branching <- simulate_controls_branching_end_node(tree, cur_node, time_after, 5)
        }
        number_branching <- c(number_branching, N_branchings)
        all_controls_branching <- c(all_controls_branching, controls_branching)
      }
    }
  }
}

df_lesion_branching = data.frame("value"=number_branching, type = "lesion")
df_control_branching = data.frame("value"=all_controls_branching, type = "control")
print(nrow(df_lesion_branching))
print(nrow(df_control_branching))

df_result_branching = rbind(df_control_branching, df_lesion_branching)
control_mean_branching = mean(df_result_branching[df_result_branching$type=="control",]$value)
lesion_mean_branching = mean(df_result_branching[df_result_branching$type=="lesion",]$value)

p_violin_branching<-(ggplot(df_result_branching, aes(y=value, x=type, fill=type)) 
                  + geom_violin() 
                  + geom_jitter(size=0.5,col="gray35")
                  + theme_bw()
                  + ylab(paste("Number of branchings in ", time_after, " mutations", sep=""))
                  + scale_y_continuous(trans="log",breaks = c(1,10,100))
                  + scale_fill_manual(values=c("grey","coral3"))
                  + theme(legend.position = "none"))

result_counts_branching = as.data.frame(table(df_result_branching))
result_counts_branching$Frequency = ifelse(result_counts_branching$type=="control", result_counts_branching$Freq/sum(result_counts_branching[result_counts_branching$type=="control",]$Freq), result_counts_branching$Freq/sum(result_counts_branching[result_counts_branching$type=="lesion",]$Freq))
p_hist_branching<-(ggplot(result_counts_branching, aes(x=value, y=Frequency, fill=type)) 
                + geom_histogram(position="dodge", stat="identity") 
                + theme_bw()
                + xlab(paste("Number of branchings in ", time_after, " mutations", sep=""))
                + scale_fill_manual(values=c("darkgrey","coral3"))
                + geom_vline(xintercept = lesion_mean_branching, col="coral3", linetype="dashed")
                + geom_vline(xintercept = control_mean_branching, col="darkgrey", linetype="dashed")
)

result_counts_branching_compact=result_counts_branching[as.numeric(result_counts_branching$value) <10,]
result_counts_branching_compact$value = as.numeric(result_counts_branching_compact$value)
result_counts_branching_compact = rbind(result_counts_branching_compact, c("10+", "control", sum(result_counts_branching[result_counts_branching$type == "control" & as.numeric(result_counts_branching$value)>=10,]$Freq),sum(result_counts_branching[result_counts_branching$type == "control" & as.numeric(result_counts_branching$value)>=10,]$Frequency)))
result_counts_branching_compact = rbind(result_counts_branching_compact, c("10+", "lesion", sum(result_counts_branching[result_counts_branching$type == "lesion" & as.numeric(result_counts_branching$value)>=10,]$Freq),sum(result_counts_branching[result_counts_branching$type == "lesion" & as.numeric(result_counts_branching$value)>=10,]$Frequency)))

result_counts_branching_compact$Frequency=as.numeric(result_counts_branching_compact$Frequency)
p_hist_branching_compact<-(ggplot(result_counts_branching_compact, aes(x=factor(value, levels=c("1","2","3","4","5","6","7","8","9","10+")), y=Frequency, fill=type)) 
                   + geom_histogram(position="dodge", stat="identity") 
                   + theme_bw()
                   + xlab(paste("Number of branchings in ", time_after, " mutations", sep=""))
                   + scale_fill_manual(values=c("darkgrey","coral3"), name="")
                   + geom_vline(xintercept = lesion_mean_branching, col="coral3", linetype="dashed")
                   + geom_vline(xintercept = control_mean_branching, col="darkgrey", linetype="dashed")
)


ks_branching = ks.test(df_result_branching[df_result_branching$type=="lesion",]$value,df_result_branching[df_result_branching$type=="control",]$value)
pval_branching = format(ks_branching$p.value,digits=3,scientific = T)

p_cumul_branching<- (ggplot(df_result_branching, aes(x=value, color=type, y = 1 - after_stat(y))) 
                  + stat_ecdf(geom="step") 
                  + theme_bw()
                  + xlab(paste("Number of branchings in ", time_after, " mutations", sep=""))
                  + scale_color_manual(values=c("darkgrey","coral3"))
                  + ylab("1-ecdf")
                  + annotate(geom="text",x=20, y=0.5, label=paste("pval=",pval_branching))
                  + theme(legend.position="none")              
)

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
tiff(filename="Phylogeny_properties_lesion_vs_control_PVVS_not_prebirth_notMPN_with_pval.jpg", width=40, height=20, res=300, units='cm')
layout_matrix = rbind(c(1, 1, 2,3), c(4,4,5,6), c(7,7,8,9))
grid.arrange(p_hist_leaves, p_cumul_leaves,p_violin_leaves,p_hist_branching, p_cumul_branching,p_violin_branching, p_hist_duration, p_cumul_duration,p_violin_duration, nrow=3,layout_matrix=layout_matrix)
dev.off()

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
png(filename="Branchngs_in_200muts_lesion_vs_control_PVVS_not_prebirth_notMPN_with_pval.png", width=20, height=8, res=300, units='cm')
grid.arrange(p_hist_branching_compact, p_cumul_branching,nrow=1)
dev.off()

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
png(filename="Distances_lesion_vs_control_PVVS_not_prebirth_notMPN_with_pval.png", width=20, height=8, res=300, units='cm')
grid.arrange(p_hist_duration, p_cumul_duration,nrow=1)
dev.off()

################## Drivers in the same node ######################################

lesion_drivers <- NULL
all_controls_drivers<-NULL
for (sample in hspc_samples){
  print(sample)
  sample_id = hspc_ids[which(hspc_samples == sample)]
  sample_drivers = drivers[drivers$Sample_ID == sample_id,]
  sample_drivers_nodes = sample_drivers$node
  print("Number of sample drivers")
  print(nrow(sample_drivers))
  tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  table_sample = str_split_i(sample, "_standard",1)
  table_sample = str_split_i(table_sample, "_noMixed",1)
  print(table_sample)
  sample_pvvs = df_muts[df_muts$Sample_ID == table_sample,]
  sample_pvvs <- sample_pvvs[sample_pvvs$Type =="PVV" & sample_pvvs$Class == "PASS",]
  print(head(sample_pvvs))
  if (nrow(sample_pvvs) != 0){
    end_nodes = sample_pvvs$lesion_repair_node
    print(end_nodes)
    start_nodes = sample_pvvs$lesion_node
    for (i in 1:length(end_nodes)){
      cur_node = start_nodes[i]
      if (nodeheight(tree, cur_node) > 60){
        lesion_drivers <- c(lesion_drivers, sum(cur_node %in% sample_drivers_nodes))
        controls <- simulate_controls_drivers(tree, cur_node, sample_drivers_nodes, 10)
        print("Simulated drivers")
        print(controls)
        all_controls_drivers <- c(all_controls_drivers, controls)
      }
      
    }
  }
}
print(sum(lesion_drivers)/length(lesion_drivers))
print(sum(all_controls_drivers)/length(all_controls_drivers))
n_drivers_lesion = sum(lesion_drivers)
n_lesion = length(lesion_drivers)
n_drivers_control = sum(all_controls_drivers)
n_controls = length(all_controls_drivers)

drivers_result = data.frame("counts" = c(n_drivers_lesion,n_lesion - n_drivers_lesion, n_drivers_control, n_controls-n_drivers_control), "Freq" = c(n_drivers_lesion/n_lesion,(n_lesion - n_drivers_lesion)/n_lesion, n_drivers_control/n_controls, (n_controls - n_drivers_control)/n_controls),"type" = c("lesion","lesion","control","control"), "driver_occurence" = c("yes","no","yes","no"))
df_for_chisq <- matrix(c(sum(lesion_drivers),sum(all_controls_drivers),length(lesion_drivers),length(all_controls_drivers)),nrow=2,ncol=2)

p_drivers_lesion_node<-(ggplot(drivers_result, aes(x=type, y=Freq, fill=driver_occurence))
  + geom_bar(stat="identity")
  + theme_bw()
  + scale_fill_manual(values = c("grey", "red"))
  + annotate(geom="text", x=1, y=1.05, label = paste(n_drivers_control, "/", n_controls))
  + annotate(geom="text", x=2, y=1.05, label = paste(n_drivers_lesion, "/", n_lesion))
  + annotate(geom="text", x=1.5, y=1.15, label = paste("pval =",round(chisq.test(df_for_chisq)$p.val, 2)))
  + geom_segment(aes(x = 1, y = 1.1, xend = 2, yend = 1.1))
  + ggtitle("Driver in the lesion start node")
)


################## Drivers any node higher ######################################
lesion_drivers_up <- NULL
all_controls_drivers_up<-NULL
for (sample in hspc_samples){
  print(sample)
  sample_id = hspc_ids[which(hspc_samples == sample)]
  sample_drivers = drivers[drivers$Sample_ID == sample_id,]
  sample_drivers_nodes = sample_drivers$node
  print("Number of sample drivers")
  print(nrow(sample_drivers))
  tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  table_sample = str_split_i(sample, "_standard",1)
  table_sample = str_split_i(table_sample, "_noMixed",1)
  print(table_sample)
  sample_pvvs = df_muts[df_muts$Sample_ID == table_sample,]
  sample_pvvs <- sample_pvvs[sample_pvvs$Type =="PVV" & sample_pvvs$Class == "PASS",]
  print(head(sample_pvvs))
  if (nrow(sample_pvvs) != 0){
    end_nodes = sample_pvvs$lesion_repair_node
    print(end_nodes)
    start_nodes = sample_pvvs$lesion_node
    for (i in 1:length(end_nodes)){
      cur_node = start_nodes[i]
      if (nodeheight(tree, cur_node) > 60){
        print(sample_drivers_nodes)
        print(getAncestors(tree, cur_node))
        print(intersect(sample_drivers_nodes, getAncestors(tree, cur_node)))
        driver_up <- ifelse(length(intersect(sample_drivers_nodes, getAncestors(tree, cur_node))) >= 1, 1,0)
        lesion_drivers_up <- c(lesion_drivers_up, driver_up)
        dist2root = length(getAncestors(tree, cur_node, type = "all"))
        controls <- simulate_controls_drivers_up(tree, cur_node, sample_drivers_nodes, dist2root, 10)
        print("Simulated drivers")
        print(controls)
        all_controls_drivers_up <- c(all_controls_drivers_up, controls)
      }
      
    }
  }
}
print(sum(lesion_drivers_up)/length(lesion_drivers_up))
print(sum(all_controls_drivers_up)/length(all_controls_drivers_up))
n_drivers_lesion_up = sum(lesion_drivers_up)
n_lesion_up = length(lesion_drivers_up)
n_drivers_control_up= sum(all_controls_drivers_up)
n_controls_up = length(all_controls_drivers_up)

drivers_result_up = data.frame("counts" = c(n_drivers_lesion_up,n_lesion_up - n_drivers_lesion_up, n_drivers_control_up, n_controls_up-n_drivers_control_up), "Freq" = c(n_drivers_lesion_up/n_lesion_up,(n_lesion_up - n_drivers_lesion_up)/n_lesion_up, n_drivers_control_up/n_controls_up, (n_controls_up - n_drivers_control_up)/n_controls_up),"type" = c("lesion","lesion","control","control"), "driver_occurence_up" = c("yes","no","yes","no"))
df_for_chisq_up <- matrix(c(sum(lesion_drivers_up),sum(all_controls_drivers_up),length(lesion_drivers_up),length(all_controls_drivers_up)),nrow=2,ncol=2)

p_drivers_lesion_node_up<-(ggplot(drivers_result_up, aes(x=type, y=Freq, fill=driver_occurence_up))
                        + geom_bar(stat="identity")
                        + theme_bw()
                        + scale_fill_manual(values = c("grey", "red"))
                        + annotate(geom="text", x=1, y=1.05, label = paste(n_drivers_control_up, "/", n_controls_up))
                        + annotate(geom="text", x=2, y=1.05, label = paste(n_drivers_lesion_up, "/", n_lesion_up))
                        + annotate(geom="text", x=1.5, y=1.15, label = paste("pval =",round(chisq.test(df_for_chisq_up)$p.val, 2)))
                        + geom_segment(aes(x = 1, y = 1.1, xend = 2, yend = 1.1))
                        + ggtitle("Driver up to the lesion start node")
)

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
tiff(filename="Driver_presence_lesion_vs_control_PVVS_not_prebirth_notMPN_with_pval.jpg", width=20, height=10, res=300, units='cm')
grid.arrange(p_drivers_lesion_node, p_drivers_lesion_node_up, nrow=1)
dev.off()

################## Branching for different times #################
regime = "start"
times = c(20,50,100,200,300,400,500,600,700,800,900,1000)
control_means<-NULL
lesion_means<-NULL
control_sizes<-NULL
lesion_sizes<-NULL
for (time_after in times){
  number_branching <- NULL
  lesion_nodes <- NULL
  lesion_repair_nodes <- NULL
  all_controls_branching<-NULL
  for (sample in hspc_samples){
    print(sample)
    tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
    tree <- di2multi(read.tree(tree_file_path))
    table_sample = str_split_i(sample, "_standard",1)
    table_sample = str_split_i(table_sample, "_noMixed",1)
    print(table_sample)
    sample_pvvs = df_muts[df_muts$Sample_ID == table_sample,]
    sample_pvvs <- sample_pvvs[sample_pvvs$Type =="PVV" & sample_pvvs$Class == "PASS",]
    print(head(sample_pvvs))
    lesion_nodes <- c(lesion_nodes, sample_pvvs$lesion_node)
    lesion_repair_nodes <- c(lesion_repair_nodes, sample_pvvs$lesion_repair_node)
    print(nrow(sample_pvvs))
    if (nrow(sample_pvvs) != 0){
      end_nodes = sample_pvvs$lesion_repair_node
      print(end_nodes)
      start_nodes = sample_pvvs$lesion_node
      for (i in 1:length(end_nodes)){
        if (regime == "start"){
          cur_node = start_nodes[i]
        }else{
          cur_node = end_nodes[i]
        }
        check_longest_offspring = FALSE
        all_cur_node_leaves = get_all_node_leaves(cur_node, tree)
        for (leaf in all_cur_node_leaves){
          if (nodeheight(tree,leaf)-nodeheight(tree,cur_node)>time_after){
            check_longest_offspring = TRUE
          }
        }
        if (check_longest_offspring){
          N_branchings = length(get_branching_points_in_time(tree, cur_node, time_after))
          print("Observed")
          print(N_branchings)
          if (regime == "start"){
            controls_branching <- simulate_controls_branching_start_node(tree, cur_node, time_after, 5)
          }else{
            controls_branching <- simulate_controls_branching_end_node(tree, cur_node, time_after, 5)
          }
          number_branching <- c(number_branching, N_branchings)
          all_controls_branching <- c(all_controls_branching, controls_branching)
        }
    }
  }
}

df_lesion_branching = data.frame("value"=number_branching, type = "lesion")
df_control_branching = data.frame("value"=all_controls_branching, type = "control")
print(nrow(df_lesion_branching))
print(nrow(df_control_branching))
df_result_branching = rbind(df_control_branching, df_lesion_branching)
control_mean_branching = mean(df_result_branching[df_result_branching$type=="control",]$value)
control_sample_size = length(all_controls_branching)
lesion_sample_size = length(number_branching)
lesion_mean_branching = mean(df_result_branching[df_result_branching$type=="lesion",]$value)
control_means <- c(control_means, control_mean_branching)
lesion_means <- c(lesion_means, lesion_mean_branching)
control_sizes <- c(control_sizes, control_sample_size)
lesion_sizes <- c(lesion_sizes, lesion_sample_size)
}
df_by_time = data.frame("Time"=c(times, times), "Branchings" = c(lesion_means, control_means), "Sample_size" = c(lesion_sizes,control_sizes),"Type"=c(rep("lesion",length(times)), rep("control",length(times))))

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
tiff(filename=paste("Branching_vs_time_lesion_", regime, "_node.jpg"), width=15, height=15, res=300, units='cm')
(ggplot(df_by_time, aes(x=Time, y=Branchings, col=Type, group=Type)) 
                   + geom_point()
                    + geom_line()
                   + theme_bw()
                   + xlab("Time in mutations")
                  + ylab("Mean number of branchings")
                  + geom_text(aes(label=Sample_size), nudge_x = -30, size=3)
                   + scale_color_manual(values=c("darkgrey","coral3"))
                  + ggtitle(paste("Lesion ",regime, " node", sep=""))
)
dev.off()
################################
hspc_samples = c("PD4781","PD5117","PD5163","PD5179","PD5182","PD5847","PD6629", "PD6646","PD7271","PD9478","PD34493_FINAL_noMixed","PD41276_v2_FINAL_rho05_noMixed","PD41305_v3_FINAL_noMixed","AX001_4_01_standard_rho01","CB001_3_01_standard_rho01","CB002_2_01_standard_rho01","KX001_4_01_standard_rho01","KX002_2_01_standard_rho01","KX003_5_01_standard_rho01","KX004_5_01_standard_rho01","KX007_2_01_standard_rho01","KX008_2_01_standard_rho01","KX009_1_01_standard_rho01","KX010_1_01_standard_rho01","SX001_5_01_standard_rho01")

distances <- NULL
for (sample in hspc_samples){
  print(sample)
  tree_file_path = paste0(data_dir,"/trees/tree_",sample,".tree")
  tree <- di2multi(read.tree(tree_file_path))
  for (node in c(Ntip(tree)+2:Nnode(tree))){
    print(node)
    node_ancestors <- getAncestors(tree, node)
    print(node_ancestors)
    if (length(node_ancestors) >= 2){
      if (length(node_ancestors) > 4){node_ancestors = node_ancestors[1:4]}
      for (ancestor in node_ancestors[1:(length(node_ancestors)-1)]){
      #for (ancestor in node_ancestors[1]){
        dist <- nodeheight(tree, node) - nodeheight(tree, ancestor)
        distances<-c(distances, dist)
      }
    }
  }
}

print((distances))
hist(distances, breaks=100)

hist.data = hist(distances, plot=F, breaks=100)
hist.data$counts=hist.data$counts+1
hist.data$counts = log10(hist.data$counts)
plot(hist.data, ylab='log(Frequency)')




#####################
node = 533
EdgeCols <- rep("black", Nedge(tree))
EdgeCols[which(tree$edge[,2] == node,)]<-"red"
plot(tree, edge.color =EdgeCols)
