library(data.table)
library(dplyr)
library(purrr)
library(stats4)
library(tidyr)
library(ggplot2)
library(ggridges)

# Inputs:
# - norm_mut: data.frame with columns chr, pos (normal mutations)
# - pvv_mut:  data.frame with columns chr, pos (PVV mutations)

win_size <- 1e7
chr_list <- paste0("chr", 1:23)

# Step 1: Estimate local mutation rate per window
estimate_local_rate <- function(norm_mut, win_size=1e7) {
  norm_mut <- norm_mut %>% mutate(bin = floor(Pos / win_size))
  rates <- norm_mut %>% group_by(Chrom, bin) %>% summarise(rate = n()) %>% ungroup()
  rates <- rates %>% group_by(Chrom) %>% mutate(rate = rate / sum(rate)) %>% ungroup()
  return(rates)
}

# Step 2: Simulate PVV inheritance and recombination
simulate_pvv <- function(rates, d, lambda, win_size=1e7) {
  sim_result <- list()
  n_chr <- length(unique(rates$Chrom))
  chr_names <- unique(rates$Chrom)
  # Each homolog starts with its own active state (both initially carry PVVs)
  chr_states <- lapply(chr_names, function(chr) {
    chr_rates <- rates %>% filter(Chrom == !!chr)
    n_bins <- nrow(chr_rates)
    list(
      chr = chr,
      bin = chr_rates$bin,
      rate = chr_rates$rate,
      mat = rep(TRUE, n_bins),
      pat = rep(TRUE, n_bins)
    )
  })
  # Simulate d divisions with recombination before each
  for (i in 1:d) {
    for (idx in seq_along(chr_states)) {
      state <- chr_states[[idx]]
      n_bins <- length(state$rate)
      
      # Recombination before division
      n_recomb <- rpois(1, lambda)
      if (n_recomb >= n_bins){
        n_recomb = n_bins - 1
      }
      if (n_recomb > 0) {
        breakpoints <- sort(sample(1:(n_bins - 1), n_recomb, replace=FALSE))
        segments <- cut(1:n_bins, breaks=c(0, breakpoints, n_bins), labels=FALSE)
        new_pat <- state$pat
        new_mat <- state$mat
        
        # Alternate exchange status across segments: first is exchanged, second is not, etc.
        exchange_flags <- rep(c(TRUE, FALSE), length.out = max(segments))
        
        for (s in unique(segments)) {
          if (exchange_flags[s]) {
            # If exchange flag is TRUE, swap this segment between homologs
            new_pat[segments == s] <- state$mat[segments == s]
            new_mat[segments == s] <- state$pat[segments == s]
          }
        }
        # Update state after recombination
        state$pat <- new_pat
        state$mat <- new_mat
      }
      
      # Random segregation: 50% chance to lose one homolog
      if (runif(1) < 0.5) {
        state$pat <- rep(FALSE, n_bins)
      } else {
        state$mat <- rep(FALSE, n_bins)
      }
      
      chr_states[[idx]] <- state
      #message("division=", d)
      #message("sum_mat=",sum(state$pat))
      #message("sum_pat",sum(state$mat))
    }
  }
  
  # Normalize global mutation rates over all chromosomes
  global_rate <- unlist(lapply(chr_states, function(s) s$rate * (s$pat | s$mat)))
  if (sum(global_rate) == 0) {
    return(bind_rows(lapply(chr_states, function(s) tibble(chr=s$chr, bin=s$bin, count=rep(0, length(s$bin))))))
  }
  global_prob <- global_rate / sum(global_rate)
  counts <- rmultinom(1, size=54, prob=global_prob)[,1]
  # Assign counts back per chromosome/bin
  k <- 1
  for (idx in seq_along(chr_states)) {
    state <- chr_states[[idx]]
    n_bins <- length(state$rate)
    chr_states[[idx]]$count <- counts[k:(k+n_bins-1)]
    k <- k + n_bins
  }
  # Combine result
  sim_result <- bind_rows(lapply(chr_states, function(s) {
    tibble(chr=s$chr, bin=s$bin, count=s$count)
  }))
  return(sim_result)
}

# Step 3: Define likelihood function
likelihood <- function(d, lambda, rates, observed_pvvs, n_sim=100) {
  sims <- lapply(seq_len(n_sim), function(i) {
    sim_result <- simulate_pvv(rates, d, lambda)
    sim_result$sim <- i  # add a simulation ID column
    return(sim_result)
  })
  print(length(sims))
  print("done")
  mean_sim_counts <- bind_rows(sims) %>% group_by(sim,chr) %>% summarise(sum_sim=sum(count)) %>% group_by(chr) %>% summarise(mean_sim=mean(sum_sim))
  
  sim_chr_counts <- bind_rows(sims) %>%
    group_by(sim, chr) %>%
    summarise(sum_sim = sum(count), .groups = "drop")
  
  # For each simulation, count how many chromosomes had at least one mutation
  n_chr_with_mut <- sim_chr_counts %>%
    filter(sum_sim > 0) %>%
    group_by(sim) %>%
    summarise(n_chr = n(), .groups = "drop")
  
  # Now compute the mean number of such chromosomes across simulations
  mean_n_chr_with_mut <- mean(n_chr_with_mut$n_chr)  
  
  
  print(head(mean_sim_counts))
  message("sum simulat",sum(mean_sim_counts$mean_sim))
  obs <- observed_pvvs %>% mutate(bin = floor(Pos / win_size)) %>% group_by(Chrom) %>% summarise(obs=n())
  colnames(obs) <- c("chr",  "obs")
  print(sum(obs$obs))
  merged <- left_join(obs, mean_sim_counts, by=c("chr")) %>% mutate(mean_sim = replace_na(mean_sim, 1e-4))
  print(merged)
  print(sum(merged$obs))
  print(sum(merged$mean_sim))
  #ll <- sum(dpois(merged$obs, lambda=merged$mean_sim, log=TRUE))
  #ll <- mean_n_chr_with_mut
  #return(ll)
  n_chr_with_mut$LAD <- d
  n_chr_with_mut$lamda <- lambda
  return(n_chr_with_mut)
}

# Step 4: Grid search or optimization
fit_parameters <- function(rates, observed_pvvs) {
  results <- list()
  
  for (d in 1:6) {
#  for (d in 1:4) {
   for (lambda in seq(0, 3, 1)) {
#    for (lambda in c(0)) {
      cat("d=", d, "lambda=", lambda,  "\n")
      ll <- likelihood(d, lambda, rates, observed_pvvs, n_sim=50)
      message("ll=",ll)
      results[[length(results) + 1]] <- tibble(d = d, lambda = lambda, ll = ll)
    }
  }
  
  results_df <- bind_rows(results)
  return(results_df)
}

# Step 5: Produce the whole distribution from simulations
full_simul_results <- function(rates, observed_pvvs) {
  full_result <- NULL
  for (d in 1:5) {
    #  for (d in 1:4) {
    for (lambda in seq(0, 3, 1)) {
      #    for (lambda in c(0)) {
      cat("d=", d, "lambda=", lambda,  "\n")
      ll <- likelihood(d, lambda, rates, observed_pvvs, n_sim=50)
      full_result <- rbind(full_result, ll)
    }
  }
  return(full_result)
}

# Example usage:
#norm_mut <- read.table("norm_mut.tsv", header=TRUE)
#pvv_mut <- read.table("pvv_mut.tsv", header=TRUE)

interesting_node = 418
df_multi = read.table("/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/trees_and_mutations/mutations_filtered.tsv", header=T, sep="\t")
pvv_mut = df_multi[df_multi$Sample_ID == "PX001_2_01" & df_multi$Type == "PVV"& df_multi$lesion_node == interesting_node,]
print(nrow(pvv_mut))


load("/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/trees_and_mutations/MF/temp/annotated_mut_set_PX001_2_01_standard_rho01")
norm_mut = filtered_muts$COMB_mats.tree.build$mat

rates <- estimate_local_rate(norm_mut)

fit <- fit_parameters(rates, pvv_mut)
print(fit)
# Find best parameters
#best_row <- fit %>% filter(ll == max(ll, na.rm = TRUE)) %>% slice(1)
#best_params <- c(d = best_row$d, lambda = best_row$lambda)

# Plot likelihoods
ggplot(fit, aes(x = lambda, y = factor(d), fill = ll)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Log-Likelihood") +
  labs(title = "Log-Likelihood Surface", x = "Lambda", y = "Divisions (d)") +
  theme_minimal()


##############################
df_multi = read.table("/home/mandrianova/Lab/Burst_kinetics/Mike/from_Mike/trees_and_mutations/mutations_filtered.tsv", header=T, sep="\t")

df_multi_PVV_chemo = df_multi[df_multi$Sample_ID == "PX001_2_01" & df_multi$Type == "PVV"& df_multi$lesion_node == interesting_node,]
df_multi_PVV_chemo$chr = str_split_i(df_multi_PVV_chemo$Chrom_pos, "-",1)
df_multi_PVV_chemo_agg = as.data.frame(df_multi_PVV_chemo %>% group_by(chr) %>% count())
colnames(df_multi_PVV_chemo_agg) <- c("Chrom", "MultA")
Obs<-sum(as.numeric(df_multi_PVV_chemo_agg$MultA>0.1))


final_result <- full_simul_results(rates, pvv_mut)

setwd("/home/mandrianova/Lab/Burst_kinetics/Mike/for_Mike/results/")
jpeg("Node418_LAD_with_recombination.jpeg", width = 13, height = 10, units = "cm", res = 300)

(ggplot(final_result[final_result$lamda!=3 & final_result$LAD != 5,], aes(x=n_chr,y=factor(LAD, levels = c(4,3,2,1)), fill=as.factor(lamda))) 
       #+ geom_density_ridges(alpha=0.6,scale=0.8,stat = "binline", bins=5))
       + geom_density_ridges(alpha=0.6)
  + geom_vline(xintercept = Obs, size=1)
  + theme(plot.title = element_text(hjust = 0.5))
  + theme_bw()
  #+ scale_fill_manual(values=colors_list[2:6])
  #+ scale_color_manual(values=colors_list[2:6])
  + ylab("LAD")
  + theme(legend.position = "right")
  + xlab("N chroms with multiallelic sites")
  + scale_fill_manual(values=c("#B42D73","#EE316B","#FFB137","#FBDFB6"), name="recombinations\n per chromosome")
)

dev.off()



