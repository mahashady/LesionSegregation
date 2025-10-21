library(ggplot2)
library(stringr)
library(lsa)
library(tidyverse)
library(pheatmap)
library(grid)

df_expected1 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/enriched.Platinum.SBS_less10.expected_triallelic.csv", header=TRUE, sep=",")
df_expected1$enriched.Platinum.SBS17_less10.expected <- df_expected1$expected_count/sum(df_expected1$expected_count)
df_expected2 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/NONenriched.Platinum.SBS_less10.expected_triallelic.csv", header=TRUE, sep=",")
df_expected2$NONenriched.Platinum.SBS17_less10.expected <- df_expected2$expected_count/sum(df_expected2$expected_count)
# df_expected3 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/enriched.Platinum.SBS_more10.expected_triallelic.csv", header=TRUE, sep=",")
# df_expected3$enriched.Platinum.SBS17_more10.expected <- df_expected3$expected_count/sum(df_expected3$expected_count)
# df_expected4 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/NONenriched.Platinum.SBS_more10.expected_triallelic.csv", header=TRUE, sep=",")
# df_expected4$NONenriched.Platinum.SBS17_more10.expected <- df_expected4$expected_count/sum(df_expected4$expected_count)
# df_expected_list <- list(df_expected1, df_expected2, df_expected3, df_expected4)
df_expected_list <- list(df_expected1, df_expected2)
# df_expected_list <- list(df_expected3, df_expected4)
df_expected <- reduce(df_expected_list, full_join, by = "context")
print(head(df_expected))

df_observed1 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/multi_spectrum/multi_spectrum_enriched.Platinum.SBS17_less10.tsv", header=TRUE, sep=",")[,c("mutation","Type","N_muts")]
df_observed1$enriched.Platinum.SBS17_less10.observed <- df_observed1$N_muts/sum(df_observed1$N_muts)
df_observed2 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/multi_spectrum/multi_spectrum_NONenriched.Platinum.SBS17_less10.tsv", header=TRUE, sep=",")[,c("mutation","N_muts")]
df_observed2$NONenriched.Platinum.SBS17_less10.observed <- df_observed2$N_muts/sum(df_observed2$N_muts)
# df_observed3 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/multi_spectrum/multi_spectrum_enriched.Platinum.SBS17_more10.tsv", header=TRUE, sep=",")[,c("mutation","Type","N_muts")]
# df_observed3$enriched.Platinum.SBS17_more10.observed <- df_observed3$N_muts/sum(df_observed3$N_muts)
# df_observed4 <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/multi_spectrum/multi_spectrum_NONenriched.Platinum.SBS17_more10.tsv", header=TRUE, sep=",")[,c("mutation","N_muts")]
# df_observed4$NONenriched.Platinum.SBS17_more10.observed <- df_observed4$N_muts/sum(df_observed4$N_muts)
# df_observed_list <- list(df_observed1, df_observed2, df_observed3, df_observed4)
df_observed_list <- list(df_observed1, df_observed2)
# df_observed_list <- list(df_observed3, df_observed4)
df_observed <- reduce(df_observed_list, full_join, by = "mutation")
print(head(df_observed))

# result <- merge(df_expected[,c("context","enriched.Platinum.SBS17_less10.expected","NONenriched.Platinum.SBS17_less10.expected","enriched.Platinum.SBS17_more10.expected", "NONenriched.Platinum.SBS17_more10.expected")], df_observed[,c("mutation","Type","enriched.Platinum.SBS17_less10.observed","NONenriched.Platinum.SBS17_less10.observed","enriched.Platinum.SBS17_more10.observed","NONenriched.Platinum.SBS17_more10.observed")], by.x="context", by.y="mutation")
result <- merge(df_expected[,c("context","enriched.Platinum.SBS17_less10.expected","NONenriched.Platinum.SBS17_less10.expected")], df_observed[,c("mutation","Type","enriched.Platinum.SBS17_less10.observed","NONenriched.Platinum.SBS17_less10.observed")], by.x="context", by.y="mutation")
# result <- merge(df_expected[,c("context","enriched.Platinum.SBS17_more10.expected","NONenriched.Platinum.SBS17_more10.expected")], df_observed[,c("mutation","Type","enriched.Platinum.SBS17_more10.observed","NONenriched.Platinum.SBS17_more10.observed")], by.x="context", by.y="mutation")
print(head(result))
df_signatures <- read.table("/workspace/projects/lesion_segregation/metastatic_tumors/results/expected_triallelic_spectrum/signatures.expected_triallelic.Platinum.csv", header=TRUE, sep=",")
df_signatures$context <- rownames(df_signatures)
print(head(df_signatures))
# result <- merge(result, df_signatures[,c("X","X25_1", "X37_1", "X31_SBS17b_0.968799_1", "X21_SBS31_0.953955_1", "X14_1")], by.x="context", by.y="X")
result <- merge(result, df_signatures[,c("X","X25_1", "X37_1", "X21_SBS31_0.953955_1", "X14_1")], by="context", by.y="X")
print(head(result))

print(colnames(result))
result <- result %>% 
  rename(
    'Carboplatin\n(E-SBS25)' = X25_1,
    'Oxaliplatin\n(E-SBS37)' = X37_1,
    'Carbo/Cisplatin\n(E-SBS21)' = 'X21_SBS31_0.953955_1',
    'Cis/Oxaliplatin\n(E-SBS14)' = X14_1,
    # 'E-SBS31' ='X31_SBS17b_0.968799_1',    
    'MAV enriched exp' = 'enriched.Platinum.SBS17_less10.expected',
    'MAV non-enriched exp' = 'NONenriched.Platinum.SBS17_less10.expected',
    'enriched samples' = 'enriched.Platinum.SBS17_less10.observed',
    'non-enriched samples' = 'NONenriched.Platinum.SBS17_less10.observed'
    )

print(head(result))

cosine_result <- NULL
for (i in 5:6){
    for (j in c(2:3,7:10)){
        cos <- cosine(result[,i],result[,j])
        cosine_result <- rbind(cosine_result,c(colnames(result)[i], colnames(result)[j], cos))
    }

}


cosine_result <- as.data.frame(cosine_result)
print((cosine_result))
colnames(cosine_result) <- c("x", "y", "cosine")
cosine_result$cosine <- as.numeric(cosine_result$cosine)
print(head(cosine_result))
# jpeg(filename="../plots/compare_multiallelic_expected_vs_observed.jpeg", width=30, height=15, res=300, units='cm')

# (ggplot(cosine_result, aes(y, x, fill= cosine)) 
#     + geom_tile()
#     + scale_fill_gradient(low="grey", high = "red")
# )
# dev.off()

# Create cosine similarity matrix
cols_x <- 5:6
cols_y <- c(2:3,7:10)

mat <- matrix(NA, nrow=length(cols_x), ncol=length(cols_y),
              dimnames=list(colnames(result)[cols_x], colnames(result)[cols_y]))

for (i in seq_along(cols_x)) {
  for (j in seq_along(cols_y)) {
    mat[i, j] <- cosine(result[, cols_x[i]], result[, cols_y[j]])
  }
}
mat <- as.data.frame(mat)
mat[["MAV\nexpected"]] <- ifelse(
  rownames(mat) == "enriched samples", mat[["MAV enriched exp"]],
  ifelse(
    rownames(mat) == "non-enriched samples", mat[["MAV non-enriched exp"]],
    NA  # better to use NA than "" for numeric data
  )
)
mat <- mat[, !names(mat) %in% c("MAV enriched exp", "MAV non-enriched exp")]
mat <- mat[, c(ncol(mat), 1:(ncol(mat)-1))]
print(mat)
# Plot heatmap with clustering
# Open JPEG device
library(pheatmap)
library(grid)

# Open JPEG device
jpeg(filename = "../plots/compare_multiallelic_expected_vs_observed_clustered.Platinum.SBS17less10.jpeg", 
     width = 20, height = 9, res = 300, units = 'cm')

# Generate heatmap as a grob
p <- pheatmap(mat,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              clustering_method = "complete",
              color = colorRampPalette(c("grey", "aquamarine4"))(50),
              display_numbers = TRUE,
              number_format = "%.2f",
              fontsize_number = 10,
              number_color = "black",
              angle_col = 45,
              cellwidth = 50,   
              cellheight = 50,
              silent = TRUE)

# Start a new blank grid page
grid.newpage()

# Create layout with margins
pushViewport(viewport(layout = grid.layout(3, 3,
    widths = unit.c(unit(1, "cm"), unit(1, "null"), unit(1, "cm")),
    heights = unit.c(unit(1, "cm"), unit(1, "null"), unit(1, "cm"))
)))

# Plot heatmap in the center of the layout (row 2, col 2)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(p$gtable)
popViewport()

# Add legend title (above colorbar)
# Approximate position – adjust x/y as needed
grid.text("cosine", 
          x = unit(0.84, "npc"), 
          y = unit(0.95, "npc"), 
          gp = gpar(fontsize = 12, fontface = "bold"))

# Close device
dev.off()