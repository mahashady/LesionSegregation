root_folder <- "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors"
df_pretreatment <- read.table(paste0(root_folder, "/data/pre_biopsy_drugs.tsv"), header=TRUE, sep="\t")
df_pretreatment_agg_name <- aggregate(df_pretreatment$name, list(df_pretreatment$patientIdentifier), FUN=paste, collapse= ",")
colnames(df_pretreatment_agg_name)<- c("patientIdentifier", "name")
df_pretreatment_agg_mechanism <- aggregate(df_pretreatment$mechanism, list(df_pretreatment$patientIdentifier), FUN=paste, collapse= ",")
colnames(df_pretreatment_agg_mechanism)<- c("patientIdentifier", "mechanism")
df_pretreatment_agg = merge(df_pretreatment_agg_name,df_pretreatment_agg_mechanism, by="patientIdentifier")

df_pretreatment_agg$Immunotherapy = ifelse(grepl("Anti-PD-1|Anti-PD-L1", df_pretreatment_agg$mechanism), 1, 0)
df_pretreatment_agg$Platinum = ifelse(grepl("Platinum", df_pretreatment_agg$mechanism), 1, 0)
df_pretreatment_agg$Alkylating = ifelse(grepl("Alkylating", df_pretreatment_agg$mechanism), 1, 0)

df_pretreatment_agg$Cisplatin = ifelse(grepl("Cisplatin", df_pretreatment_agg$name), 1, 0)
df_pretreatment_agg$Carboplatin = ifelse(grepl("Carboplatin", df_pretreatment_agg$name), 1, 0)
df_pretreatment_agg$Oxaliplatin = ifelse(grepl("Oxaliplatin", df_pretreatment_agg$name), 1, 0)
df_pretreatment_agg$Capecitabine = ifelse(grepl("Capecitabine", df_pretreatment_agg$name), 1, 0)

print(nrow(df_pretreatment_agg))
write.csv(df_pretreatment_agg, file = paste0(root_folder, "/results/pre_biopsy_drugs.agg.txt"), row.names=F)
