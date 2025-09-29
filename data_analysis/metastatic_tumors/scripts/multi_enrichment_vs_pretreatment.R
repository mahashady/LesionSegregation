library(dplyr)
library(ggplot2)

df_pretreatment = read.table("../results/pre_biopsy_drugs.agg.txt", header = T, sep=",")
#head(df_pretreatment)

df_multi = read.table("../results/ALL_Hartwig_bi_multi.txt", header=TRUE, sep="\t")
df_multi$patientIdentifier = substr(as.character(df_multi$sample), start= 1,  stop= nchar(as.character(df_multi$sample))-1)
#print(head(df_multi))

multi = merge(df_multi, df_pretreatment, by = "patientIdentifier")
print("Merged all")
print(nrow(multi))
#multi$n_multi_exp = (multi$n_biallelic)/3000000000*(multi$n_biallelic)
multi$n_multi_exp = (multi$n_biallelic)**2/(3*10**9)
multi$enrichment = multi$n_multi/multi$n_multi_exp
multi$enriched = ifelse(multi$enrichment > 10, "Yes", "No")

multi_selected = multi[(multi$Immunotherapy+multi$Platinum+multi$Alkylating)==1,]
print("N samples immuno")
print(nrow(multi_selected[multi_selected$Immunotherapy == 1, ]))
print("N samples platinum")
print(nrow(multi_selected[multi_selected$Platinum == 1, ]))
print("N samples alkyl")
print(nrow(multi_selected[multi_selected$Alkylating == 1, ]))
multi_selected$therapy = ifelse(multi_selected$Immunotherapy == 1, "Anti-PD-1|Anti-PD-L1", ifelse(multi_selected$Platinum == 1, "Platinum", "Alkylating"))
multi_selected$therapy_type = ifelse(multi_selected$therapy == "Anti-PD-1|Anti-PD-L1", "Immunotherapy", "Chemotherapy")
print(head(multi_selected))
multi_nonenriched = multi_selected[multi_selected$enriched == "No",]
multi_enriched = multi_selected[multi_selected$enriched == "Yes",]
write.table(multi_enriched, file = "../results/all_enriched.chemo_alkyl_immuno.txt", row.names=F, quote=F, sep="\t")
write.table(multi_nonenriched, file = "../results/all_NONenriched.chemo_alkyl_immuno.txt", row.names=F, quote=F, sep="\t")

multi_agg = as.data.frame(multi_selected %>% group_by(therapy_type, therapy, enriched) %>% summarise(count=n()) %>% mutate(perc=count/sum(count)))
multi_agg_enriched = multi_agg[multi_agg$enriched=="Yes",]
print(multi_agg)

jpeg(filename="../plots/Hartwig_multiallelic_enriched_vs_treatment.jpeg", width=10, height=7, res=300, units='cm')
(ggplot(multi_agg_enriched, aes(x=therapy, y=perc, fill=factor(therapy, levels = c("Alkylating", "Platinum", "Anti-PD-1|Anti-PD-1"))))
  + geom_bar(stat="identity")
  + facet_grid(~therapy_type,scales="free", space='free')
  + theme_bw()
  + scale_fill_manual(name="", values=c("cornflowerblue","darkorchid4","cornsilk3"))
  + xlab("")
  + theme(legend.position = "none")
  + geom_text(aes(label=count), vjust=-0.2)
  + ylab("Enriched samples, %")
  + ylim(0,0.15)
  )
dev.off()