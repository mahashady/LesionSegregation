library("betareg")
library("dplyr")
library(flexmix)
library(stringr)

### This script runs fitting of vaf distribution by 2 beta distributions only for samples MRCA=1 and all 4 cells survived (HMM_multi state2 > 95%)
### We DONT want to exclude multi-allelic sites as this can bias our estimations. So we take all the variants and calculate vaf as 1-ref/sum(a+t+g+c)

args = commandArgs(trailingOnly=TRUE)
driver = args[1]
print(driver)
map_nucl=c("A"="T", "T"="A","G"="C","C"="G")
map_strand=c("+"="-", "-"="+")

fit_betamix <- function(df){
#    if (nrow(df) > 5000){
#  	  df_subset = sample_n(df, 5000, replace=F)
#    }else{
#	    df_subset = df}
    df_subset=df
    result = c(NA,NA,NA,NA)  
    tryCatch(
    {model <- betamix(vaf ~ 1 | 1, data = df_subset, k = 2, nstart=3)
    print("Model is ready")
    print(length(coef(model)))
    if (length(coef(model)) == 4){
      mu1 <- plogis(coef(model)[1])
      mu2 <- plogis(coef(model)[2])
      phi1 <- exp(coef(model)[3])
      phi2 <- exp(coef(model)[4])
      a1 <- mu1 * phi1
      b1 <- (1 - mu1) * phi1
      a2 <- mu2 * phi2
      b2 <- (1 - mu2) * phi2
      p1 <- prior(model$flexmix)[1]
      p2 <- prior(model$flexmix)[2]  
        if (mu1 < mu2){
          weight_subclon = p1
          weight_clon = p2
          mu_subclon= mu1
          mu_clon = mu2
        }else{
          weight_subclon = p2
          weight_clon = p1
          mu_subclon = mu2
          mu_clon = mu1}
    }else{
      mu1 <- plogis(coef(model)[1])
      phi1 <- exp(coef(model)[3])
      a1 <- mu1 * phi1
      b1 <- (1 - mu1) * phi1
      p1 <- prior(model$flexmix)[1]
      weight_clon = p1
      weight_subclon = 0
      mu_clon = mu1
      mu_subclon=NA
    }
    print("HERE!")
    result = c(mu_subclon, mu_clon, weight_subclon, weight_clon)
    print(result)
    }, error = function(e){
            message('Caught an error!')
        }
  )
  return(result)
}

df_MRCA <- read.table("../../MRCA/results/Summary_divisions_with_symmetrical_no_tetra.txt", header=TRUE, sep=",")
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
  file = paste("../../data/mutations/", sample, ".nodMat", sep="")
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
    result = fit_betamix(df_cur)
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
write.table(result_df, file = paste("fit_beta_mixture_1div_4cells_full.txt", sep=""), quote = F, sep = "\t",
            row.names = F, col.names = T)