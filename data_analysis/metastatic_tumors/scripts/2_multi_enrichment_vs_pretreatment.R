library(dplyr)
library(ggplot2)
library(grid)

df_pretreatment <- read.table("../results/pre_biopsy_drugs.agg.txt", header = TRUE, sep = ",")
df_multi <- read.table("../results/ALL_Hartwig_bi_multi.txt", header = TRUE, sep = "\t")

df_multi$patientIdentifier <- substr(
  as.character(df_multi$sample),
  start = 1,
  stop = nchar(as.character(df_multi$sample)) - 1
)

multi <- merge(df_multi, df_pretreatment, by = "patientIdentifier")

print("Merged all")
print(nrow(multi))

multi$n_multi_exp <- (multi$n_biallelic)^2 / (3 * 10^9)
multi$enrichment <- multi$n_multi / multi$n_multi_exp
multi$enriched <- ifelse(multi$enrichment > 10, "Yes", "No")

multi_selected <- multi[
  (multi$Immunotherapy + multi$Platinum + multi$Alkylating) == 1,
]

print("N samples immuno")
print(nrow(multi_selected[multi_selected$Immunotherapy == 1, ]))

print("N samples platinum")
print(nrow(multi_selected[multi_selected$Platinum == 1, ]))

print("N samples alkyl")
print(nrow(multi_selected[multi_selected$Alkylating == 1, ]))

multi_selected$therapy <- ifelse(
  multi_selected$Immunotherapy == 1,
  "Anti-PD-1|Anti-PD-L1",
  ifelse(multi_selected$Platinum == 1, "Platinum", "Alkylating")
)

multi_selected$therapy_type <- ifelse(
  multi_selected$therapy == "Anti-PD-1|Anti-PD-L1",
  "Immunotherapy",
  "Chemotherapy"
)

multi_nonenriched <- multi_selected[multi_selected$enriched == "No", ]
multi_enriched <- multi_selected[multi_selected$enriched == "Yes", ]

write.table(
  multi_enriched,
  file = "../results/enrichment/all_enriched.chemo_alkyl_immuno.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

write.table(
  multi_nonenriched,
  file = "../results/enrichment/all_NONenriched.chemo_alkyl_immuno.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

multi_agg <- as.data.frame(
  multi_selected %>%
    group_by(therapy_type, therapy, enriched) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(perc = count / sum(count)) %>%
    ungroup()
)

multi_agg_enriched <- multi_agg[multi_agg$enriched == "Yes", ]

# Fisher tests vs immunotherapy
subset_alkyl <- multi_selected[
  multi_selected$therapy %in% c("Alkylating", "Anti-PD-1|Anti-PD-L1"),
]

subset_alkyl$therapy <- factor(
  subset_alkyl$therapy,
  levels = c("Alkylating", "Anti-PD-1|Anti-PD-L1")
)

subset_alkyl$enriched <- factor(
  subset_alkyl$enriched,
  levels = c("Yes", "No")
)

subset_platinum <- multi_selected[
  multi_selected$therapy %in% c("Platinum", "Anti-PD-1|Anti-PD-L1"),
]

subset_platinum$therapy <- factor(
  subset_platinum$therapy,
  levels = c("Platinum", "Anti-PD-1|Anti-PD-L1")
)

subset_platinum$enriched <- factor(
  subset_platinum$enriched,
  levels = c("Yes", "No")
)

print(table(subset_alkyl$therapy, subset_alkyl$enriched))
fisher_alkyl <- fisher.test(
  table(subset_alkyl$therapy, subset_alkyl$enriched),
  alternative = "greater"
)

print(table(subset_platinum$therapy, subset_platinum$enriched))
fisher_platinum <- fisher.test(
  table(subset_platinum$therapy, subset_platinum$enriched),
  alternative = "greater"
)

label_alkyl <- paste0("p = ", signif(fisher_alkyl$p.value, 2))
label_platinum <- paste0("p = ", signif(fisher_platinum$p.value, 2))

therapy_levels <- c("Alkylating", "Platinum", "Anti-PD-1|Anti-PD-L1")

multi_agg_enriched$therapy <- factor(
  multi_agg_enriched$therapy,
  levels = therapy_levels
)

# P-value bracket annotation dataframe
annot_df <- data.frame(
  xmin = c("Alkylating", "Platinum"),
  xmax = c("Anti-PD-1|Anti-PD-L1", "Anti-PD-1|Anti-PD-L1"),
  y = c(0.18, 0.205),
  label = c(label_alkyl, label_platinum),
  stringsAsFactors = FALSE
)

annot_df$xmin <- factor(annot_df$xmin, levels = therapy_levels)
annot_df$xmax <- factor(annot_df$xmax, levels = therapy_levels)

# Boxed group labels under x-axis
group_df <- data.frame(
  x = c(1.5, 3),
  y = c(-0.035, -0.035),
  label = c("Chemotherapy", "Immunotherapy")
)

jpeg(
  filename = "../plots/Hartwig_multiallelic_enriched_vs_treatment.jpeg",
  width = 10,
  height = 7,
  res = 300,
  units = "cm"
)

ggplot(
  multi_agg_enriched,
  aes(x = therapy, y = perc, fill = therapy)
) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_bw() +
  scale_fill_manual(
    name = "",
    values = c(
      "Alkylating" = "cornflowerblue",
      "Platinum" = "darkorchid4",
      "Anti-PD-1|Anti-PD-L1" = "cornsilk3"
    )
  ) +
  xlab("") +
  ylab("Fraction of samples") +
  theme(
    legend.position = "none",
    #axis.text.x = element_text(angle = 30, hjust = 1),
    plot.margin = margin(5.5, 5.5, 40, 5.5)
  ) +
  geom_text(
    aes(label = count),
    vjust = -0.4,
    size = 3
  ) +

  # p-value brackets
  geom_segment(
    data = annot_df,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = annot_df,
    aes(x = xmin, xend = xmin, y = y - 0.007, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = annot_df,
    aes(x = xmax, xend = xmax, y = y - 0.007, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = annot_df,
    aes(x = xmin, y = y + 0.008, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    hjust = -0.6
  ) +

  # boxed therapy type labels below the axis
  geom_label(
    data = group_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    label.size = 0.25,
    label.padding = unit(0.15, "lines"),
    fill = "white"
  ) +
  coord_cartesian(ylim = c(-0.055, 0.22), clip = "off")

dev.off()