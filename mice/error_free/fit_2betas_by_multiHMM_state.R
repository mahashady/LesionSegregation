library("betareg")
library("dplyr")
library(flexmix)
library(stringr)

### For samples with more than 1 multiHMM state, we want to fit mixture of betas separately for each state
### This will include 1st division (if cells died) or divisions later than 1
###We DONT want to exclude multi-allelic sites as this can bias our estimations. So we take all the variants and calculate vaf as (1-ref)/sum(a+t+g+c)


fit_betamix <- function(df){
    if (nrow(df) > 5000){
  	  df_subset = sample_n(df, 5000, replace=F)
    }else{
	    df_subset = df}
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

df_samples_C3H = read.table("../../data/c3h.tumourCollateInfo.tab", sep=",", header=T)
df_samples_CAST = read.table("../../data/cast.tumourCollateInfo.tab", sep=",", header=T)
df_samples = rbind(df_samples_C3H, df_samples_CAST)
all_samples_list = df_samples$nodId


df_MRCA <- read.table("../../MRCA/results/Summary_divisions_with_symmetrical_no_tetra.txt", header=TRUE, sep=",")
asym_samples = df_MRCA[df_MRCA$division != 0,]$sample
print(length(asym_samples))

n_A0_list <- NULL
prop_A0_list <- NULL
w_subclon_A0_list <- NULL
w_clon_A0_list <- NULL
mu_subclon_A0_list <- NULL
mu_clon_A0_list <- NULL
n_A1_list <- NULL
prop_A1_list <- NULL
w_subclon_A1_list <- NULL
w_clon_A1_list <- NULL
mu_subclon_A1_list <- NULL
mu_clon_A1_list <- NULL
n_A2_list <- NULL
prop_A2_list <- NULL
w_subclon_A2_list <- NULL
w_clon_A2_list <- NULL
mu_subclon_A2_list <- NULL
mu_clon_A2_list <- NULL
i= 0



for (sample in asym_samples){
    i = i + 1
    print(i)
    print(sample)
    file = paste("../../data/mutations/", sample, ".nodMat", sep="")
    file_multiHMM = paste("../../MRCA/results/HMM/output_HMM/", sample, ".hmm.PloidyHMM", sep="")
    df_multiHMM = read.table(file_multiHMM)
    print(nrow(df_multiHMM))
    file_inputHMM = paste("../../MRCA/results/HMM/input_HMM/", sample, ".hmm", sep="")
    df_inputHMM = read.table(file_inputHMM)
    df_inputHMM$multiHMM_state = df_multiHMM$V1
    df_inputHMM$ref = str_split_i(df_inputHMM$V3, "_", 1)
    df_inputHMM$alt = str_split_i(df_inputHMM$V4, "_", 2)
    df_inputHMM$muId = paste(df_inputHMM$V1, ":",df_inputHMM$V2, "_", df_inputHMM$ref, "/", df_inputHMM$alt, sep="")

    df_cur = read.table(file, sep=",", header=T)
    print(nrow(df_cur))
    df_cur = df_cur[df_cur$isIndel!= 1,] #select only SNPs
    print(nrow(df_cur))
  
    df_cur = merge(df_cur, df_inputHMM[,c("muId", "multiHMM_state")], by="muId")
    df_cur = df_cur[df_cur$chr!= 'X',] #select autosomes]
    df_cur$vaf = 1-df_cur$refCount/(df_cur$au+df_cur$cu+df_cur$tu+df_cur$gu)


    df_state1 = df_cur[df_cur$multiHMM_state == "A0",]
    df_state2 = df_cur[df_cur$multiHMM_state == "A1",]
    df_state3 = df_cur[df_cur$multiHMM_state == "A2",]

    n_A0_list = c(n_A0_list, nrow(df_state1))
    n_A1_list = c(n_A1_list, nrow(df_state2))
    n_A2_list = c(n_A2_list, nrow(df_state3))

    prop_A0_list = c(prop_A0_list, (nrow(df_state1)/nrow(df_cur)))
    prop_A1_list = c(prop_A1_list, (nrow(df_state2)/nrow(df_cur)))
    prop_A2_list = c(prop_A2_list, (nrow(df_state3)/nrow(df_cur)))


    if ((nrow(df_state1)/nrow(df_cur)) < 0.05){
	    mu_subclon_A0_list = c(mu_subclon_A0_list, NA)
	    mu_clon_A0_list = c(mu_clon_A0_list, NA)
	    w_subclon_A0_list = c(w_subclon_A0_list, NA)
      w_clon_A0_list = c(w_clon_A0_list, NA)
    }else{
      print("Modeling for state A0")
      result_A0 = fit_betamix(df_state1)
      print(result_A0)
    	mu_subclon_A0_list = c(mu_subclon_A0_list, result_A0[1])
	    mu_clon_A0_list = c(mu_clon_A0_list, result_A0[2])
    	w_subclon_A0_list = c(w_subclon_A0_list, result_A0[3])
      w_clon_A0_list = c(w_clon_A0_list, result_A0[4])
      print(result_A0[3])
    }
  
    if ((nrow(df_state2)/nrow(df_cur)) < 0.05){
      mu_subclon_A1_list = c(mu_subclon_A1_list, NA)
      mu_clon_A1_list = c(mu_clon_A1_list, NA)
	    w_subclon_A1_list = c(w_subclon_A1_list, NA)
  	  w_clon_A1_list = c(w_clon_A1_list, NA)
    }else{
      print("Modeling for state A1")
      result_A1 = fit_betamix(df_state2)
      print(result_A1)
	    mu_subclon_A1_list = c(mu_subclon_A1_list, result_A1[1])
	    mu_clon_A1_list = c(mu_clon_A1_list, result_A1[2])
	    w_subclon_A1_list = c(w_subclon_A1_list, result_A1[3])
        w_clon_A1_list = c(w_clon_A1_list, result_A1[4])
        print(result_A1[3])
    }

    if ((nrow(df_state3)/nrow(df_cur)) < 0.05){
	    mu_subclon_A2_list = c(mu_subclon_A2_list, NA)
	    mu_clon_A2_list = c(mu_clon_A2_list, NA)
	    w_subclon_A2_list = c(w_subclon_A2_list, NA)
  	  w_clon_A2_list = c(w_clon_A2_list, NA)
    }else{
      print("Modeling for state A2")
      result_A2 = fit_betamix(df_state3)
      print(result_A2)
	    mu_subclon_A2_list = c(mu_subclon_A2_list, result_A2[1])
	    mu_clon_A2_list = c(mu_clon_A2_list, result_A2[2])
	    w_subclon_A2_list = c(w_subclon_A2_list, result_A2[3])
      w_clon_A2_list = c(w_clon_A2_list, result_A2[4])
      print(result_A2[3])
    }

}  
print(length(n_A0_list))
print(length(prop_A0_list))
print(length(mu_clon_A0_list))
print(length(mu_subclon_A0_list))
print(length(w_subclon_A0_list))
print(length(w_clon_A0_list))
print(length(n_A1_list))
print(length(prop_A1_list))
print(length(mu_clon_A1_list))
print(length(mu_subclon_A1_list))
result_df = data.frame(sample = asym_samples, 
    n_A0 = n_A0_list, prop_A0 = prop_A0_list, mu_subclon_A0 = mu_subclon_A0_list, mu_clon_A0 = mu_clon_A0_list, w_subclon_A0 = w_subclon_A0_list, w_clon_A0 = w_clon_A0_list,
    n_A1 = n_A1_list, prop_A1 = prop_A1_list, mu_subclon_A1 = mu_subclon_A1_list, mu_clon_A1 = mu_clon_A1_list, w_subclon_A1 = w_subclon_A1_list, w_clon_A1 = w_clon_A1_list,
    n_A2 = n_A2_list, prop_A2 = prop_A2_list, mu_subclon_A2 = mu_subclon_A2_list, mu_clon_A2 = mu_clon_A2_list, w_subclon_A2 = w_subclon_A2_list, w_clon_A2 = w_clon_A2_list
    )
head(result_df)


setwd("../../error_free/results/")
write.table(result_df, file = "fit_beta_mixture_whole_genome_by_multiHMM_state.txt", quote = F, sep = "\t",
            row.names = F, col.names = T)