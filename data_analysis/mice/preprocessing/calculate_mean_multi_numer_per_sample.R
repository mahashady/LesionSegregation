
df <- read.table("../data/c3h.tumourCollateInfo.tab", header=TRUE, sep=",")
print(head(df))
samples_C3H <- df$nodId
print(samples_C3H)

n_multi_all <- NULL
for (sample in samples_C3H){
    file_name <- paste0("../LAD/results/HMM/input_HMM/", sample, ".hmm")
    df_sample <- read.table(file_name, sep=" ")
    print(head(df_sample))
    n_multi <- nrow(df_sample[df_sample$V8 == "M",])
    n_multi_all <- c(n_multi_all, n_multi)
}

print(n_multi_all)
print(min(n_multi_all))
print(max(n_multi_all))
print(median(n_multi_all))