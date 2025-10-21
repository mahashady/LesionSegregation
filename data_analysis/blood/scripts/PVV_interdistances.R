library(dplyr)
library(ggplot2)
library(stringr)
library(GenomeInfoDb)
library(tibble)

# Toggle: include edge distances or not
include_edges <- FALSE

# PVV data
df <- read.table("/workspace/projects/lesion_segregation/blood/data/mutations_filtered.tsv", header=TRUE, sep="\t")
print(head(df))
print(nrow(df))
df <- df[df$Sample_ID == "PX001_2_01",]
print(nrow(df))
df <- df[(df$Type == "PVV" & df$Class == "PASS")|(df$Type == "MAV" & df$Class == "PASS"),]
print(nrow(df))
df <- df[df$lesion_node == 418,]
print(nrow(df))
df$chrom = paste0("chr",str_split_i(df$Chrom_pos, "-", 1))
df$pos = as.numeric(str_split_i(df$Chrom_pos, "-", 2))
df <- df[,c("chrom", "pos", "Type", "Class")]
print(nrow(df))

# Chromosome lengths
chrom_info <- getChromInfoFromUCSC("hg38")

chrom_lengths <- chrom_info %>%
  filter(!grepl("_", chrom)) %>%              # remove alternative contigs
  filter(chrom %in% paste0("chr", c(1:22, "X", "Y"))) %>%
  select(chrom, size) %>%
  deframe()  # Convert to named vector
print(chrom_lengths)


# Function to compute normalized distances for one chromosome
get_norm_distances <- function(df_chr, chrom_len, include_edges) {
  df_chr <- df_chr %>% arrange(pos)
  n <- nrow(df_chr)
  if (n < 2) return(numeric(0))

  positions <- df_chr$pos
  print(positions)
  internal <- diff(positions)

  if (include_edges) {
    from_start <- positions[1]
    to_end <- chrom_len - positions[n]
    all_distances <- c(from_start, internal, to_end)
  } else {
    all_distances <- internal
  }
  print(all_distances)
  all_distances / (chrom_len * (1/(n+1)))
}

# Compute for all chromosomes
# df_chr1 <- df[df$chrom == "chr1",]
# test <- get_norm_distances(df_chr1, chrom_lengths["chr1"], include_edges)
# print(test)

normalized_distances <- df %>%
  filter(chrom %in% names(chrom_lengths)) %>%
  group_by(chrom) %>%
  group_map(~ get_norm_distances(.x, chrom_lengths[.y$chrom[1]], include_edges)) %>% 
  unlist()

print(normalized_distances)
df_dist <- data.frame(norm_distance = normalized_distances)
print(head(df_dist))

# Save plot
jpeg(
  filename = paste0("/workspace/projects/lesion_segregation/blood/plots/PVV_interdistances", 
       if (include_edges) "_with_edges" else "", 
       ".jpeg"),
  width = 10, height = 7, res = 300, units = 'in'  # 'in' is usually better with high DPI
)

ggplot(df_dist, aes(x = norm_distance)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(
    x = paste("Normalized Distance", if (include_edges) "(with edges)" else ""),
    y = "Frequency",
    title = paste("Distribution of Normalized PVV Distances", if (include_edges) "(with edges)" else "")
  ) +
  theme_minimal()

dev.off()