library(mixtools)
library("dplyr")
library(flexmix)
library(stringr)

### This script runs fitting of vaf distribution by 2 beta distributions only for samples MRCA=1 and all 4 cells survived (HMM_multi state2 > 95%)
### Before fitting we subset mutations corresponing to the driver of interest by expression_quantile, gene strand and substitution type.
### We DONT want to exclude multi-allelic sites as this can bias our estimations. So we take all the variants and calculate vaf as 1-ref/sum(a+t+g+c)

args = commandArgs(trailingOnly=TRUE)
i = as.integer(args)
drivers <- c("7:145859242_T/A", "7:145859242_T/C", "6:37548568_A/T", "11:14185624_T/A")
driver = drivers[i]
print(driver)
map_nucl=c("A"="T", "T"="A","G"="C","C"="G")
map_strand=c("+"="-", "-"="+")

fit_norm_mixture <- function(df){
  model <- normalmixEM(df$vaf, k = 2) #fit vaf distribution as a mixture of two normal distributions
  w1 <- sum(model$posterior[,1]>0.5)/nrow(model$posterior)
  w2 <- sum(model$posterior[,2]>0.5)/nrow(model$posterior)
  m1 <- model$mu[1]
  m2 <- model$mu[2]
  print(m1)
  print(m2)
  if (is.na(m2)){
    mu_subclon <- NA
    mu_clon <- m1
    w_subclon <- NA
    w_clon <- w2
  }else{
    if (m1 < m2){
      mu_subclon <- m1
      mu_clon <- m2
      w_subclon <- w1
      w_clon <- w2
    }else{
      mu_subclon <- m2
      mu_clon <- m1
      w_subclon <- w2
      w_clon <- w1
    }

  }
  result = c(mu_subclon, mu_clon, w_subclon, w_clon)
  return(result)
}

df_samples_C3H = read.table("../../data/c3h.tumourCollateInfo.tab", sep=",", header=T)
df_samples_CAST = read.table("../../data/cast.tumourCollateInfo.tab", sep=",", header=T)
df_samples = rbind(df_samples_C3H, df_samples_CAST)
all_samples_list = df_samples$nodId


df_MRCA <- read.table("../../MRCA/results/Summary_divisions_with_symmetrical_no_tetra.txt", header=TRUE, sep=",")
samples_div1_cells4 <- df_MRCA[df_MRCA$division == 1 & df_MRCA$prop_multi_hmm_state3 > 0.95 & df_MRCA$mice_line == "C3H",]$sample
print("Number of samples=")
print(length(samples_div1_cells4))

df_drivers_annotated <- read.table("/workspace/projects/lesion_segregation/mice/data/Drivers_annotated.csv", header=TRUE, sep=",")
driver_expr_q <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$expression_q
driver_strand <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$Gene_strand
driver_chr <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$chr
driver_pos <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$pos
driver_ref <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$refN
driver_alt <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$altN
driver_gene <- df_drivers_annotated[df_drivers_annotated$mutID == driver,]$Gene_name


w_subclon_list <- NULL
w_clon_list <- NULL
mu_subclon_list <- NULL
mu_clon_list <- NULL
n_mut_list <- NULL

for (sample in samples_div1_cells4){
  print(sample)
  file = paste("../../data/mutations_vs_genes/", sample, ".with_gene_annot.nodMat", sep="")
  df_cur = read.table(file, sep=",", header=T)
  print("all mutations")
  print(nrow(df_cur))
  df_cur = df_cur[df_cur$isIndel!= 1,] #select only SNPs
  df_cur = df_cur[df_cur$chr!= 'X',] #select autosomes
  print("SNPs and autosomes")
  print(nrow(df_cur))
  df_cur$vaf = 1-df_cur$refCount/(df_cur$au+df_cur$cu+df_cur$tu+df_cur$gu)
  df_cur <- df_cur[df_cur$expression_q == driver_expr_q,] # subset mutations with the same expression quantile as driver
  print("corresponding expression")
  print(nrow(df_cur))
  df_cur <- df_cur[(df_cur$geneStrand == driver_strand & df_cur$refN == driver_ref & df_cur$altN == driver_alt)|(df_cur$geneStrand == map_strand[driver_strand] & df_cur$refN == map_nucl[driver_ref] & df_cur$altN == map_nucl[driver_alt]),]
  print("final set")
  print(nrow(df_cur))
  n_mut_list <- c(n_mut_list, nrow(df_cur))
  #print(head(df_cur))
  if (nrow(df_cur) > 100){
    print("Fitting")
    result = fit_norm_mixture(df_cur)
    mu_subclon_list = c(mu_subclon_list, result[1])
	  mu_clon_list = c(mu_clon_list, result[2])
    w_subclon_list = c(w_subclon_list, result[3])
    w_clon_list = c(w_clon_list, result[4])
  }else{
    mu_subclon_list = c(mu_subclon_list, NA)
	  mu_clon_list = c(mu_clon_list, NA)
    w_subclon_list = c(w_subclon_list, NA)
    w_clon_list = c(w_clon_list, NA)

  }
}  
result_df = data.frame(sample = samples_div1_cells4, n_mut = n_mut_list,
    mu_subclon = mu_subclon_list, mu_clon = mu_clon_list, w_subclon = w_subclon_list, w_clon = w_clon_list)
head(result_df)


setwd("../../error_free/results/")
write.table(result_df, file = paste("fit_normal_mixture_1div_4cells_", driver_gene, "_", driver_chr, "_", driver_pos, "_", driver_ref, "_", driver_alt,".txt", sep=""), quote = F, sep = "\t",row.names = F, col.names = T)