library(ggplot2)
library(dplyr)
library(ggside)

divisions_all <- read.table("../results/Summary_LAD_with_symmetrical_no_mixtures.txt", header=TRUE, sep = ",")
divisions_not_sym <- divisions_all[divisions_all$LAD_ML_segm != "0",]
divisions <- divisions_all[divisions_all$LAD_ML_segm != "0" & divisions_all$LAD_ML_segm != "1",]
plot_state1_2<-(ggplot(divisions, aes(x=segments_HMM_multi_state_M1_prop, y = segments_HMM_multi_state_M2_prop,col=LAD_ML_segm)) +
                geom_segment(aes(x = 1/4, y = 0, xend = 1/4, yend = 2/4), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 0, y = 2/4, xend = 1/4, yend = 2/4), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 9/16, y = 0, xend = 9/16, yend = 6/16), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 0, y = 6/16, xend = 9/16, yend = 6/16), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 49/64, y = 0, xend = 49/64, yend = 14/64), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 0, y = 14/64, xend = 49/64, yend = 14/64), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 225/256, y = 0, xend = 225/256, yend = 30/256), linetype="dashed",col="darkgrey") +
                geom_segment(aes(x = 0, y = 30/256, xend = 225/256, yend = 30/256), linetype="dashed",col="darkgrey") +
                geom_point(x=1/4, y=2/4, col="grey40") +
                geom_point(x=9/16, y=6/16, col="grey40") +
                geom_point(x=49/64, y=14/64, col="grey40") +
                geom_point(x=225/256, y=30/256, col="grey40") +
                geom_point(size=1) +
                annotate("text", x=0.16, y=0.48, label="II division") +
                annotate("text", x=0.46, y=0.35, label="III division") +
                annotate("text", x=0.69, y=0.2, label="IV division") +
                annotate("text", x=0.81, y=0.1, label="V division") +
                theme_bw()  +
                xlab("proportion of segments in M1 state of HMM_multi") +
                ylab("proportion of segments in M2 state of HMM_multi") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                scale_x_continuous(expand = c(0, 0.01), breaks = c(0, 1/4,9/16,49/64,225/256,1.2), labels = c("0", parse(text="(1/2)^2"), parse(text="(3/4)^2"), parse(text="(7/8)^2"), parse(text="(15/16)^2"),"")) +
                scale_y_continuous(expand = c(0, 0.01), breaks = c(0, 2/4,6/16,14/64,30/256,1.2), labels = c("0", parse(text="2*(1/2)*(1/2)"), parse(text="2*(1/4)*(3/4)"), parse(text="2*(1/8)*(7/8)"), parse(text="2*(1/16)*(15/16)"), ""))
)
ggsave(plot_state1_2, file="../plots/multi_hmm_state1_2_scatter.jpeg", width = 7, height = 7, dpi = 300)

plot_state1_distr <- (ggplot(divisions, aes(x=segments_HMM_multi_state_M1_prop)) 
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
ggsave(plot_state1_distr, file="../plots/multi_hmm_state1_distribution_C3H.jpeg", width = 9, height = 5, dpi = 300)

lad_by_line <- read.table("../results/LAD_ML_segm_by_mice_line_agg.txt", header=TRUE, sep=",")
lad_by_line <- lad_by_line %>%
  group_by(mice_line) %>%
  mutate(prop_one_driver = one_driver / sum(one_driver)) %>%
  ungroup()

MRCA_distr_by_line <- (ggplot(lad_by_line, aes(x = LAD_ML_segm, y=prop_one_driver,fill=mice_line))
                       + geom_bar(stat = "identity", position = "dodge")
                       + theme_bw()
                       + ylab("proportion of samples")
                       + scale_fill_discrete(name="mice line")
                       + xlab("LAD"))
ggsave(MRCA_distr_by_line, file="../plots/LAD_distribution_by_line.jpeg", width = 6, height = 4, dpi = 300)

divisions_not_sym$prop_multiallelic = divisions_not_sym$segments_HMM_multi_state_M2_prop + divisions_not_sym$segments_HMM_multi_state_M3_prop
divisions_not_sym$driver_id = paste(divisions_not_sym$drivers, divisions_not_sym$driver_positions, sep="_")
divisions_not_sym <- divisions_not_sym[divisions_not_sym$driver_id %in% c("Egfr_11:14185624_T/A","Braf_6:37548568_A/T", "Hras_7:145859242_T/A", "Hras_7:145859242_T/C"),]
print(head(divisions_not_sym))
driver_names <- data.frame(
  driver_id = c("Egfr_11:14185624_T/A", "Braf_6:37548568_A/T", "Hras_7:145859242_T/A", "Hras_7:145859242_T/C"),
  driver_name = c("Egfr F254I", "Braf V637E", "Hras Q61L", "Hras Q61R")
)

divisions_not_sym <- divisions_not_sym %>%
  left_join(driver_names, by = "driver_id")
  
print(head(divisions_not_sym))

df_side <- divisions_not_sym %>%
    count(driver_name, LAD_ML_segm) %>%
    group_by(driver_name) %>%
    mutate(prop = n / sum(n))

plot_multi_vs_LAD <- ggplot(
  divisions_not_sym,
  aes(x = factor(LAD_ML_segm), y = prop_multiallelic)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.4,
    col = "darkgrey"
  ) +
  geom_jitter(
    aes(col = driver_name),
    width = 0.17,
    size = 1,
    alpha = 1
  ) +
  ggside::geom_xsidecol(
    data = df_side,
      aes(
        x = factor(LAD_ML_segm),
        y = prop,
        fill = driver_name
      ),
    inherit.aes = FALSE,
    width = 0.7
  ) +
  facet_wrap(~driver_name, nrow = 2) +
  labs(
    x = "infered LAD",
    y = "Proportion of multiallelic sements",
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, colour = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "none",
    ggside.panel.scale = 0.25    
  ) + 
  scale_color_manual(values=c("Egfr F254I"="#9BCD9B", "Braf V637E"="#EE6363","Hras Q61L"="#9AC0CD", "Hras Q61R"="#9AC0CD")) +
  scale_fill_manual(values=c("Egfr F254I"="#9BCD9B", "Braf V637E"="#EE6363","Hras Q61L"="#9AC0CD", "Hras Q61R"="#9AC0CD")) +
  scale_xsidey_continuous(breaks = c(0.5), limits = c(0, 1))

ggsave(
  filename = "../plots/multi_vs_LAD.C3H.by_gene.jpeg",
  plot = plot_multi_vs_LAD,
  width = 6,
  height = 6,
  dpi = 300
)
