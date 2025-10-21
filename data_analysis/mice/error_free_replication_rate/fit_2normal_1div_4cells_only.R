library(mixtools)
library("dplyr")
library(flexmix)
library(stringr)

### This script runs fitting of vaf distribution by 2 normal distributions only for samples MRCA=1 and all 4 cells survived (HMM_multi state2 > 95%)
### We DONT want to exclude multi-allelic sites as this can bias our estimations. So we take all the variants and calculate vaf as 1-ref/sum(a+t+g+c)

args = commandArgs(trailingOnly=TRUE)
driver = args[1]
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

df_MRCA <- read.table("../LAD/results/Summary_divisions_with_symmetrical_no_mixtures.txt", header=TRUE, sep=",")
samples_div1_cells4 <- df_MRCA[df_MRCA$division == 1 & df_MRCA$prop_multi_hmm_state3 > 0.95 & df_MRCA$mice_line == "C3H",]$sample
print("Number of samples=")
print(length(samples_div1_cells4))

w_subclon_list <- NULL
w_clon_list <- NULL
mu_subclon_list <- NULL
mu_clon_list <- NULL
n_mut_list <- NULL

for (sample in samples_div1_cells4){
  print(sample)
  file = paste("../data/mutations/", sample, ".nodMat", sep="")
  df_cur = read.table(file, sep=",", header=T)
  print("all mutations")
  print(nrow(df_cur))
  df_cur = df_cur[df_cur$isIndel!= 1,] #select only SNPs
  df_cur = df_cur[df_cur$chr!= 'X',] #select autosomes
  print("SNPs and autosomes")
  print(nrow(df_cur))
  df_cur$vaf = 1-df_cur$refCount/(df_cur$au+df_cur$cu+df_cur$tu+df_cur$gu)
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


setwd("../results/")
write.table(result_df, file = paste("fit_norm_mixture_1div_4cells.txt", sep=""), quote = F, sep = "\t",
            row.names = F, col.names = T)