library(hash)
library(ggplot2)
library(jsonlite)

subset <- "Alkylating"
# subset <- "Platinum"

sig_exposures = read.table("../data/Pan_full/Pan_full.exposures.tsv", header=T)
t_sig_exposures <- as.data.frame(t(sig_exposures))
print(t_sig_exposures[1:5,1:5])
print(nrow(t_sig_exposures))
t_sig_exposures$SBS31_prop = t_sig_exposures[["21_SBS31_0.953955_1"]]/rowSums(t_sig_exposures)
t_sig_exposures$SBS17b_prop = t_sig_exposures[["19_SBS17b_0.932022_1"]]/rowSums(t_sig_exposures)
t_sig_exposures$sample <- rownames(t_sig_exposures)


# Load the JSON file with genome tribases counts into a list
json_data <- fromJSON("../data/genome_counts_tribases.json")

# If it's a named vector (like a 96-mutation matrix), convert to data frame
genome_composition <- data.frame(
  context = names(json_data),
  rate = unlist(json_data)
)

colnames(genome_composition) <- c("context", "counts")
print(head(genome_composition))

muttype_level <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T" )
possible_mutations <- paste(substr(muttype_level,1,1),substr(muttype_level,3,3), substr(muttype_level,7,7), ">",substr(muttype_level,5,5), sep="")


df_enriched = read.table("../results/all_enriched.chemo_alkyl_immuno.txt", sep=",", header=TRUE)
print(head(df_enriched))
df_non_enriched = read.table("../results/all_NONenriched.chemo_alkyl_immuno.txt", sep=",", header=TRUE)
df_non_enriched = df_non_enriched[df_non_enriched$n_multi > 0,] #select NONenriched that have at least one multi-allelic sites to create the spectrum
print(nrow(df_enriched))
print(nrow(df_non_enriched))

df_enriched <- merge(df_enriched, t_sig_exposures, by="sample")
df_non_enriched <- merge(df_non_enriched, t_sig_exposures, by="sample")
# enriched_list = df_enriched[df_enriched[[subset]] == 1 & df_enriched$SBS17b_prop < 0.1,]$patientIdentifier
# non_enriched_list = df_non_enriched[df_non_enriched[[subset]] == 1 & df_non_enriched$SBS17b_prop < 0.1,]$patientIdentifier
enriched_list = df_enriched[df_enriched[[subset]] == 1,]$patientIdentifier
non_enriched_list = df_non_enriched[df_non_enriched[[subset]] == 1,]$patientIdentifier

json_samples_list <- list("enriched" = enriched_list,"NONenriched" = non_enriched_list)
# write_json(json_samples_list, paste0("/workspace/projects/lesion_segregation/metastatic_tumors/results/samples_lists.", subset, ".SBS_less10.json"), pretty = TRUE)
write_json(json_samples_list, paste0("../results/samples_lists.", subset, ".json"), pretty = TRUE)

print(length(enriched_list))
print(length(non_enriched_list))


nucls <- c("T", "C")
h <- hash()
h[["T"]] <- "A"
h[["A"]] <- "T"
h[["G"]] <- "C"
h[["C"]] <- "G"

reverse_compl <- function(context) {
  reverse_compl_context = paste(h[[substr(context,3,3)]],h[[substr(context,2,2)]],h[[substr(context,1,1)]], sep="")
  return(reverse_compl_context)
}


create_mut_1 <- function(row) {
  context <- as.character(row[5])
  alt_nucleotide <- substr(row[6], 1, 1)
  if (substr(context, 2, 2) %in% nucls) {
    mut <- paste(context, ">", alt_nucleotide, sep = "")
  } else {
    rev_compl <- reverse_compl(context)
    mut <- paste(rev_compl, ">", h[[alt_nucleotide]], sep = "")
  }
  return(mut)
}

create_mut_2 <-function(row) {
  context <- as.character(row[5])
  alt_nucleotide <- substr(row[7], 1, 1)
  if (substr(context, 2, 2) %in% nucls) {
    mut <- paste(context, ">", alt_nucleotide, sep = "")
  } else {
    rev_compl <- reverse_compl(context)
    mut <- paste(rev_compl, ">", h[[alt_nucleotide]], sep = "")
  }
  return(mut)
}

mut_dict2spectrum <- function(mut_dict, possible_mutations) {
  mut_df <- as.data.frame(table(mut_dict),stringsAsFactors=FALSE)
  colnames(mut_df) <- c("mutation", "N_muts")
  if (nrow(mut_df) != 96) {
    absent_mutations <- setdiff(possible_mutations, mut_df$mutation)
    for (abs_cont in absent_mutations){
      mut_df <- rbind(mut_df, c(abs_cont, 0))
      }
  }
  mut_df$context = substr(mut_df$mutation, 1,3)
  mut_df$nucleotide = paste(substr(mut_df$mutation, 2,2),">", substr(mut_df$mutation, 5,5), sep="")
  mut_df$Type = paste(substr(mut_df$context,1,1), "[",mut_df$nucleotide, "]",substr(mut_df$context,3,3), sep="")
  return(mut_df)
}


iter_by_samples <- function(list_samples_of_interest, subset){
  all_mutations <- NULL
  for (sample in list_samples_of_interest){
    file_name <- paste("../results/multi_sites_by_sample_", subset, ".chemo_alkyl_immuno/", sample, "_multi_sites.txt", sep="")
    all_sample_mutations <- NULL
    if (file.exists(file_name)){
      print(sample)
      df_sample = read.table(file_name)
      if (ncol(df_sample) == 8){
        df_sample$mut1 <- apply(df_sample, 1, create_mut_1)
        df_sample$mut2 <- apply(df_sample, 1, create_mut_2)
        #print(head(df_sample))
        all_mutations <- c(all_mutations, df_sample$mut1)
        all_mutations <- c(all_mutations, df_sample$mut2)
        all_sample_mutations <- c(all_sample_mutations, df_sample$mut1)
        all_sample_mutations <- c(all_sample_mutations, df_sample$mut2)
      }
      df_sample_mutations <- mut_dict2spectrum(all_sample_mutations, possible_mutations)
      df_sample_mutations <- merge(df_sample_mutations, genome_composition, by = "context")
      df_sample_mutations$Mutrate <- as.numeric(df_sample_mutations$N_muts)/as.numeric(df_sample_mutations$counts)
    }
  }
  return(all_mutations)
} 



all_mutations_enriched <- iter_by_samples(enriched_list, "enriched")
result_enriched <- mut_dict2spectrum(all_mutations_enriched, possible_mutations)
result_enriched$subset = "subset_enriched"
print(tail(result_enriched))
print("Enriched done")
print(paste("Number of enriched samples = ", length(enriched_list), sep = ""))
print(paste("Number of NONenriched samples = ", length(non_enriched_list), sep=""))

print(non_enriched_list)
all_mutations_non_enriched <- iter_by_samples(non_enriched_list, "NONenriched")
result_non_enriched <- mut_dict2spectrum(all_mutations_non_enriched, possible_mutations)
result_non_enriched$subset = "subset_non_enriched"
print("NON enriched done")

result = rbind(result_enriched, result_non_enriched)
result$N_muts <- as.numeric(result$N_muts)
print(nrow(result))
#print((result))


x_labels = c("A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T" )

p1 <- (ggplot(result)
  + geom_bar(aes(x=factor(Type,levels = muttype_level), y=N_muts, fill=nucleotide),stat="identity")
  + theme_bw()
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "top",legend.box = "horizontal",strip.text = element_text(size = 10))
  + guides(fill = guide_legend(nrow = 1))
  + scale_x_discrete(labels = x_labels)
  + scale_fill_manual(values = c("deepskyblue2","black","red","grey","darkolivegreen3","lightpink"))
  + theme(legend.key.size = unit(0.4, "cm"))
  + theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=6),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + xlab("")
  + facet_wrap(~subset, nrow = 2)
  + ggtitle(subset)
  )
ggsave(paste("../plots/multiallelic_spectrum_enriched_vs_NONenriched_", subset, ".jpeg", sep=""), p1,width=30, height=20, dpi=300, units='cm')

# write.csv(result_enriched, file = paste0("/workspace/projects/lesion_segregation/metastatic_tumors/results/multi_spectrum/multi_spectrum_enriched.", subset, ".SBS17_less10.tsv"), row.names=FALSE)
# write.csv(result_non_enriched, file = paste0("/workspace/projects/lesion_segregation/metastatic_tumors/results/multi_spectrum/multi_spectrum_NONenriched.", subset, ".SBS17_less10.tsv"), row.names=FALSE)
write.csv(result_enriched, file = paste0("../results/multi_spectrum/multi_spectrum_enriched.", subset, ".tsv"), row.names=FALSE)
write.csv(result_non_enriched, file = paste0("../results/multi_spectrum/multi_spectrum_NONenriched.", subset, ".tsv"), row.names=FALSE)
