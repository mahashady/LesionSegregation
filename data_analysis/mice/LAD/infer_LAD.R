library("readxl")
library(reshape2)
library(stringr)
library(viridis)
library(dplyr)

colors_list=c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
df_hmm_summary = read.table("../results/HMM/Summary.tsv", sep=" ")
colnames(df_hmm_summary) <- c("sample_id", "A_N2T_N_state5_high_vaf", "A_N2T_N_state5_low_vaf", "Median_vaf_autosomes", "Median_vaf_X", "n_sites_as_hmm_state1", "n_sites_as_hmm_state2","n_sites_as_hmm_state3","n_sites_as_hmm_state4","n_sites_as_hmm_state5", "n_sites_multi_hmm_state1", "n_sites_multi_hmm_state2", "n_sites_multi_hmm_state3", "emission_multi_hmm_state1", "emission_multi_hmm_state2", "emission_multi_hmm_state3", "emission_as_hmm_state1", "emission_as_hmm_state2", "emission_as_hmm_state3", "emission_as_hmm_state4", "emission_as_hmm_state5")
df_hmm_summary$sample <- str_split_i(df_hmm_summary$sample_id, "[.]", 1)
head(df_hmm_summary)

df_ploidy_summary = read.table("../results/HMM/Summary_subclone.tsv")
colnames(df_ploidy_summary) <- c("sample", "flexmix_fit", "ploidy")
head(df_ploidy_summary)

#c3h_df <- read_excel("../data/41586_2020_2435_MOESM2_ESM.xlsx", sheet = "C3H_tumours")
c3h_df <- read.table("/workspace/projects/lesion_segregation/mice/data/C3H_summary_with_sign_driver_positions.txt", header=TRUE, sep=",")
c3h_df$driver_positions_c3h_mutId <- c3h_df$driver_positions
colnames(c3h_df)[1] <- "sample"
c3h_df$mice_line = "C3H"
colnames(c3h_df)
ncol(c3h_df)
#cast_df <- read_excel("../data/41586_2020_2435_MOESM2_ESM.xlsx", sheet = "CAST_tumours")
cast_df <- read.table("/workspace/projects/lesion_segregation/mice/data/CAST_summary_with_sign_driver_positions.txt", header=TRUE, sep=",")
colnames(cast_df)
colnames(cast_df)[1] <- "sample"
cast_df$mice_line = "CAST"
ncol(cast_df)
df_mice_lines_summary = rbind(c3h_df, cast_df)
head(df_mice_lines_summary)

full_table = merge(df_hmm_summary, df_ploidy_summary, by = "sample")
full_table = merge(full_table, df_mice_lines_summary, by = "sample")
print(nrow(full_table))

full_table$prop_multi_hmm_state1 = full_table$n_sites_multi_hmm_state1/(full_table$n_sites_multi_hmm_state1 + full_table$n_sites_multi_hmm_state2 + full_table$n_sites_multi_hmm_state3)
full_table$prop_multi_hmm_state2 = full_table$n_sites_multi_hmm_state2/(full_table$n_sites_multi_hmm_state1 + full_table$n_sites_multi_hmm_state2 + full_table$n_sites_multi_hmm_state3)
full_table$prop_multi_hmm_state3 = full_table$n_sites_multi_hmm_state3/(full_table$n_sites_multi_hmm_state1 + full_table$n_sites_multi_hmm_state2 + full_table$n_sites_multi_hmm_state3)

full_table$prop_as_hmm_state1 = full_table$n_sites_as_hmm_state1/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state2 = full_table$n_sites_as_hmm_state2/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state3 = full_table$n_sites_as_hmm_state3/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state4 = full_table$n_sites_as_hmm_state4/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)
full_table$prop_as_hmm_state5 = full_table$n_sites_as_hmm_state5/(full_table$n_sites_as_hmm_state1 + full_table$n_sites_as_hmm_state2 + full_table$n_sites_as_hmm_state3 + full_table$n_sites_as_hmm_state4 + full_table$n_sites_as_hmm_state5)

full_table$Ploidy_by_HMM_as<-rep('Dip',length(full_table[,1])) #assign all samples as diploids
full_table$Ploidy_by_HMM_as[(apply(full_table[,c("prop_as_hmm_state1", "prop_as_hmm_state2", "prop_as_hmm_state3", "prop_as_hmm_state4", "prop_as_hmm_state5")] > 0.05, 1, sum)>= 4)]<-'Tetra' #assign tumors that have 4 or more as-HMM states (>0.05 state proportion) as tetraploid

full_table$tetraploids <- 0
full_table[(full_table$ploidy == "tetra" & full_table$flexmix_fit == "Good")|full_table$Ploidy_by_HMM_as == "Tetra",]$tetraploids <- 1


full_table$Symmetry <- rep('Asymmetric',length(full_table[,1]))
full_table$Symmetry[(apply(full_table[,c("prop_as_hmm_state1", "prop_as_hmm_state5")] < 0.05, 1,sum) >= 2) | full_table$emission_as_hmm_state1 < 0.8 | full_table$emission_as_hmm_state5 > 0.2 | full_table$ploidy == "Symmetric"]<- "Symmetric"
colnames(full_table)
print(full_table[full_table$Symmetry == "Symmetric",])
table(full_table[,c("Symmetry", "mice_line")])

#divisions = full_table[full_table$tetraploids != 1 & full_table$Symmetry != "Symmetric",]
tetraploids = full_table[full_table$tetraploids == 1 ,]
write.csv(tetraploids[,c("sample", "mice_line","prop_as_hmm_state1", "prop_as_hmm_state2", "prop_as_hmm_state3", "prop_as_hmm_state4", "prop_as_hmm_state5", "Ploidy_by_HMM_as", "ploidy", "drivers", "knownDrivers", "driver_positions", "driver_positions_c3h_mutId")], file = "../results/Summary_tetraploids.txt" , row.names = F, quote=FALSE)

#filtered_table <- full_table %>% filter(grepl(";", knownDrivers))
#table(filtered_table$tetraploids)
#table(full_table$tetraploids)

divisions = full_table[full_table$tetraploids != 1 ,]
divisions$division <- "X"
divisions[divisions$prop_multi_hmm_state1 > 0.95,]$division <- "late"
divisions[divisions$prop_multi_hmm_state1 < 0.05,]$division <- "1"
divisions[divisions$Symmetry == "Symmetric",]$division <- "0"
print(nrow(divisions))
print(colnames(divisions))

#print(table(divisions[divisions$division == "1",c("mice_line")]))
#print(table(divisions[divisions$division == "1" & divisions$prop_multi_hmm_state2 == 0,]$mice_line))
#print(table(divisions[divisions$division == "1" & divisions$prop_multi_hmm_state2 == 1,]$mice_line))
#print((divisions[divisions$Symmetry != "Symmetric" & divisions$prop_multi_hmm_state2 > 0.7 & divisions$prop_multi_hmm_state3 == 0,]))

#print(head(divisions[divisions$ploidy == "1_clone",]))
#print(nrow(divisions[divisions$ploidy == "1_clone",]))
#print(table(full_table[full_table$tetraploids == 1,]$mice_line))
#print(table(full_table$mice_line))
###division

divisions$dist_div2 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(1/4))^2+(as.numeric(x["prop_multi_hmm_state2"])-(2/4))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/4))^2)})
divisions$dist_div3 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(9/16))^2+(as.numeric(x["prop_multi_hmm_state2"])-(6/16))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/16))^2)})
divisions$dist_div4 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(49/64))^2+(as.numeric(x["prop_multi_hmm_state2"])-(14/64))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/64))^2)})
divisions$dist_div5 = apply(divisions, 1, function(x) {sqrt((as.numeric(x["prop_multi_hmm_state1"])-(225/256))^2+(as.numeric(x["prop_multi_hmm_state2"])-(30/256))^2+(as.numeric(x["prop_multi_hmm_state3"])-(1/256))^2)})
divisions$min_dist_cluser = apply(divisions, 1, function(x) {as.character(which.min(c(as.numeric(x["dist_div2"]),as.numeric(x["dist_div3"]),as.numeric(x["dist_div4"]),as.numeric(x["dist_div5"])))+1)})
divisions$division = ifelse(divisions$division == "X", divisions$min_dist_cluser, divisions$division)

print(nrow(divisions))
write.csv(divisions[,c("sample", "mice_line","knownDrivers", "drivers", "driver_positions", "driver_positions_c3h_mutId","division", "prop_as_hmm_state1", "prop_as_hmm_state2", "prop_as_hmm_state3", "prop_as_hmm_state4", "prop_as_hmm_state5", "emission_as_hmm_state1", "emission_as_hmm_state2", "emission_as_hmm_state3", "emission_as_hmm_state4", "emission_as_hmm_state5","prop_multi_hmm_state1","prop_multi_hmm_state2", "prop_multi_hmm_state3", "emission_multi_hmm_state1", "emission_multi_hmm_state2", "emission_multi_hmm_state3", "dist_div2", "dist_div3", "dist_div4", "dist_div5", "min_dist_cluser","numSnvs", "numIndels" )], file = "../results/Summary_mixed_excluded.txt" , row.names = F)

divisions_one_driver <-  divisions %>% filter(!grepl(";", drivers)) %>% filter(!drivers == "")
divisions_one_knownDriver <-  divisions %>% filter(!grepl(";", knownDrivers)) %>% filter(!knownDrivers == "")
divisions_by_line_all <- as.data.frame(table(divisions[,c("mice_line", "division")]))
colnames(divisions_by_line_all)[3] <- "all_samples"

divisions_by_line_one_driver <- as.data.frame(table(divisions_one_driver[,c("mice_line", "division")]))
colnames(divisions_by_line_one_driver)[3] <- "one_driver"
divisions_by_line_one_knownDriver <- as.data.frame(table(divisions_one_knownDriver[,c("mice_line", "division")]))
#colnames(divisions_by_line_one_knownDriver)[3] <- "one_strong_driver"
divisions_by_line <- merge(divisions_by_line_all,divisions_by_line_one_driver,by=c("mice_line", "division"))
#divisions_by_line <- merge(divisions_by_line,divisions_by_line_one_knownDriver,by=c("mice_line", "division"))
print(divisions_by_line)
write.csv(divisions_by_line, file = "../results/2024-11-11_LAD_by_mice_line_agg.txt" , row.names = F, quote=FALSE)


divisions_one_driver$driver_with_poisiton <- paste(divisions_one_driver$drivers, divisions_one_driver$driver_positions_c3h_mutId, sep="_")
divisions_by_line_by_driver_one_driver <- as.data.frame(table(divisions_one_driver[,c("mice_line", "division", "driver_with_poisiton")]))
write.csv(divisions_by_line_by_driver_one_driver, file = "../results/2024-11-11_LAD_by_mice_line_by_driver_agg.txt" , row.names = F, quote=FALSE)



divisions_many_drivers <- divisions %>% filter(grepl(";", drivers))
divisions_many_drivers$comments <- ""
divisions_many_drivers[divisions_many_drivers$division == "1" & divisions_many_drivers$prop_multi_hmm_state3 == 1,]$comments <- "4_cells_survived"
divisions_many_drivers[divisions_many_drivers$division == "1" & divisions_many_drivers$prop_multi_hmm_state2 > 0.05 & divisions_many_drivers$prop_multi_hmm_state3 > 0.05,]$comments <- "3_cells_survived"
divisions_many_drivers[divisions_many_drivers$division == "1" & divisions_many_drivers$prop_multi_hmm_state2 == 1,]$comments <- "not_all_survived+different_clone_size???"
divisions_many_drivers[divisions_many_drivers$division != "0" & divisions_many_drivers$division != "1" & divisions_many_drivers$prop_multi_hmm_state3 == 0 &  divisions_many_drivers$prop_multi_hmm_state1 > 0.05,]$comments <- "LAD=1 or LAD=2not_all_survived+different_clone_size??"
divisions_many_drivers <- divisions_many_drivers[,c("sample", "mice_line", "drivers","division", "comments")]
tetra_many_drivers <- tetraploids %>% filter(grepl(";", drivers))
tetra_many_drivers$comments <- ifelse(tetra_many_drivers$Ploidy_by_HMM_as == "Tetra", "as_by_HMM_as", "as_by_clon-subclon_HMM_as")
tetra_many_drivers$division <- "Mixture"
divisions_many_drivers <- divisions_many_drivers[,c("sample", "mice_line", "drivers","division", "comments")]
tetra_many_drivers <- tetra_many_drivers[,c("sample", "mice_line", "drivers","division", "comments")]
df_many_drivers <- rbind(divisions_many_drivers, tetra_many_drivers)
write.csv(df_many_drivers, file = "../results/2024-11-11_Summary_many_drivers.txt" , row.names = F, quote=FALSE)
