library(ggplot2)
library(stringr)
library(lsa)
library(tidyverse)
library(pheatmap)
library(grid)


df_expected1 <- read.table("../results/expected_triallelic_spectrum/enriched.Alkylating.expected_triallelic.csv", header=TRUE, sep=",")
df_expected1$enriched.Alkylating.expected <- df_expected1$expected_count/sum(df_expected1$expected_count)
df_expected2 <- read.table("../results/expected_triallelic_spectrum/NONenriched.Alkylating.expected_triallelic.csv", header=TRUE, sep=",")
df_expected2$NONenriched.Alkylating.expected <- df_expected2$expected_count/sum(df_expected2$expected_count)
df_expected_list <- list(df_expected1, df_expected2)
df_expected <- reduce(df_expected_list, full_join, by = "context")
print(head(df_expected))

df_observed1 <- read.table("../results/multi_spectrum/multi_spectrum_enriched.Alkylating.tsv", header=TRUE, sep=",")[,c("mutation","Type","N_muts")]
df_observed1$enriched.Alkylating.observed <- df_observed1$N_muts/sum(df_observed1$N_muts)
df_observed2 <- read.table("../results/multi_spectrum/multi_spectrum_NONenriched.Alkylating.tsv", header=TRUE, sep=",")[,c("mutation","N_muts")]
df_observed2$NONenriched.Alkylating.observed <- df_observed2$N_muts/sum(df_observed2$N_muts)
df_observed_list <- list(df_observed1, df_observed2)
df_observed <- reduce(df_observed_list, full_join, by = "mutation")
print(head(df_observed))

result <- merge(df_expected[,c("context","enriched.Alkylating.expected","NONenriched.Alkylating.expected")], df_observed[,c("mutation","Type","enriched.Alkylating.observed","NONenriched.Alkylating.observed")], by.x="context", by.y="mutation")
print(head(result))
df_signatures <- read.table("../results/expected_triallelic_spectrum/signatures.expected_triallelic.Alkylating.csv", header=TRUE, sep=",")
colnames(df_signatures)[1] <-  "context"
print(head(df_signatures))

result <- merge(result, df_signatures[,c("context","cyclophosphamide_557117b73fe2", "X38_SBS2_0.996907_1", "X20_SBS13_0.948838_1")], by="context")
print(head(result))


result <- result %>% 
  rename(
    cyclophosphamide = cyclophosphamide_557117b73fe2,
    SBS2 = 'X38_SBS2_0.996907_1',
    SBS13 = 'X20_SBS13_0.948838_1',
    'MAV enriched\nexpected' = 'enriched.Alkylating.expected',
    'MAV non-enriched\nexpected' = 'NONenriched.Alkylating.expected',
    'MAV\nenriched samples' = 'enriched.Alkylating.observed',
    'MAV\nnon-enriched samples' = 'NONenriched.Alkylating.observed'    
    )





print(head(result))

cosine_result <- NULL
for (i in 5:6){
    for (j in c(2:3, 7:9)){
        cos <- cosine(result[,i],result[,j])
        cosine_result <- rbind(cosine_result,c(colnames(result)[i], colnames(result)[j], cos))
    }

}

cosine_result <- as.data.frame(cosine_result)
print(head(cosine_result))
colnames(cosine_result) <- c("x", "y", "cosine")
cosine_result$cosine <- as.numeric(cosine_result$cosine)
jpeg(filename="../plots/compare_multiallelic_expected_vs_observed.Alkylating.jpeg", width=30, height=15, res=300, units='cm')

(ggplot(cosine_result, aes(y, x, fill= cosine)) 
    + geom_tile()
    + scale_fill_gradient(low="grey", high = "red")
)
dev.off()

# Create cosine similarity matrix
cols_x <- 5:6
cols_y <- c(2:3, 7:9)

mat <- matrix(NA, nrow=length(cols_x), ncol=length(cols_y),
              dimnames=list(colnames(result)[cols_x], colnames(result)[cols_y]))

for (i in seq_along(cols_x)) {
  for (j in seq_along(cols_y)) {
    mat[i, j] <- cosine(result[, cols_x[i]], result[, cols_y[j]])
  }
}

# Plot heatmap with clustering
# Open JPEG device
library(pheatmap)
library(grid)

# Open JPEG device
jpeg(filename = "../plots/compare_multiallelic_expected_vs_observed_clustered.Alkylating.jpeg", 
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