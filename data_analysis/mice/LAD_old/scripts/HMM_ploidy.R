library(HMM)
#library(flexmix)
library(mixtools)

##############################
## 1. Parse arguments & setup
##############################
args <- commandArgs(trailingOnly = TRUE)
input_fn <- args[1]
print(input_fn)
input_path <- paste("../results/input_HMM/", input_fn, sep = "")
output_dir <- "../results/HMM/HMM_ploidy/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

##############################
## 2. Load and preprocess data
##############################

data <- read.table(file=input_path, header=TRUE, sep=",") 
data$chr <- as.character(data$chr)
# select bi-allelic autosomal T>N/A>N mutations
message("Numer of all mutations = ", nrow(data))
SM <- data[data$code_multi == 'B' & (data$code_as == 'A>N' | data$code_as == 'T>N') & data$chr != 'X',]
SM$chr <- as.numeric(SM$chr)
SM$altCount <- as.numeric(SM$altCount)

# Calculate vaf
SM$vaf <- SM$altCount/(SM$altCount+SM$refCount)
#Numeric coordinate for ordering
SM$coord <- as.numeric(SM$chr) * 1e9 + SM$pos
message("Numer of input mutations = ", nrow(SM))

#####################################
## 3. Mixture model on VAF distribution
#####################################

fit <- normalmixEM(SM$vaf, k = 2) #fit vaf distribution as a mixture of two normal distributions
print(head(fit$mu))
# Count confidently assigned variants
if (fit$mu[1] < fit$mu[2]){
    n_subclonal <- sum(fit$posterior[, 1] > 0.8)
    n_clonal <- sum(fit$posterior[, 2] > 0.8)
} else {
    n_subclonal <- sum(fit$posterior[, 2] > 0.8)
    n_clonal <- sum(fit$posterior[, 1] > 0.8)
}

#####################################
## 4. Define and run HMM
#####################################

TM <- matrix(0.002, nrow = 3, ncol = 3)
diag(TM) <- 0.192

# Initial HMM with 3 possible states of asymmetry
hmm=initHMM(c("A1","A2","A3" ),c("T>N","A>N"),
            transProbs=TM, 
            emissionProbs=rbind(c(.9,.1),
                                c(.5,.5),
                                c(.1,.9))) 


# Subset clonal and subclonal mutations
if (fit$mu[1] < fit$mu[2]){
    obs_clonal <- SM[fit$posterior[, 2] > 0.8,]
    obs_subclonal <- SM[fit$posterior[, 1] > 0.8,]
} else {
    obs_clonal <- SM[fit$posterior[, 1] > 0.8,]
    obs_subclonal <- SM[fit$posterior[, 2] > 0.8,]
}
message("N clonal mutations = ", nrow(obs_clonal))
message("N subclonal mutations = ", nrow(obs_subclonal))
#run Baum-Welch and Viterbi for clonal mutations
bw_clonal <- baumWelch(hmm,obs_clonal$code_as, maxIterations = 15)
print("Clonal HMM done")
print(bw_clonal$hmm$emissionProbs)
obs_clonal$states_clonal <- viterbi(bw_clonal$hmm, obs_clonal$code_as)
print(head(obs_clonal))
#run Baum-Welch and Viterbi for subclonal mutations
bw_subclonal <- baumWelch(hmm, obs_subclonal$code_as, maxIterations = 15)
print("Subclonal HMM done")
print(bw_subclonal$hmm$emissionProbs)
obs_subclonal$states_subclonal <- viterbi(bw_subclonal$hmm,obs_subclonal$code_as)
print(head(obs_subclonal))
#####################################
## 5. Summarize HMM outputs
#####################################

#create matrix with emissions and number of sites in each state separately for clonal and subclonal mutations and write to the output
EMP <- c(
  "Emissions_bw1_bw2",
  bw_subclonal$hmm$emissionProbs[, 1],
  bw_clonal$hmm$emissionProbs[, 1]
)


STP <- c(
  "Statecounts_bw1_bw2",
  sum(obs_subclonal$states_subclonal == "A1"),
  sum(obs_subclonal$states_subclonal == "A2"),
  sum(obs_subclonal$states_subclonal == "A3"),
  sum(obs_clonal$states_clonal == "A1"),
  sum(obs_clonal$states_clonal == "A2"),
  sum(obs_clonal$states_clonal == "A3")
)

Bw_matrix <- rbind(EMP,STP)

# output HMM summary
write.table(
  Bw_matrix,
  file = file.path(output_dir, paste0(input_fn, ".BW")),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

#####################################
## 6. VAF distribution summaries
#####################################

# output boundary vafs for clonal and subclonal mutations, 
# median vaf of clonal and subclonal mutations, number of clonal and subclonal mutations

MAX <- min(
  max(SM$vaf[fit$posterior[, 1] > 0.8]),
  max(SM$vaf[fit$posterior[, 2] > 0.8])
)
MIN <- max(
  min(SM$vaf[fit$posterior[, 1] > 0.8]),
  min(SM$vaf[fit$posterior[, 2] > 0.8])
)


if (fit$mu[1] < fit$mu[2]){
    median_subclonal <- median(SM$vaf[fit$posterior[, 1] > 0.5])
    median_clonal <- median(SM$vaf[fit$posterior[, 2] > 0.5])
} else {
    median_subclonal <- median(SM$vaf[fit$posterior[, 2] > 0.5])
    median_clonal <- median(SM$vaf[fit$posterior[, 1] > 0.5])
}

Distr_properties <- c(
  "vafMed_vafboundaries_Distr_size",
  median_subclonal,
  median_clonal,
  MAX,
  MIN,
  n_subclonal,
  n_clonal
)
write.table(
  Distr_properties,
  file = file.path(output_dir, paste0(input_fn, ".Clonesize")),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

#####################################
## 7. Asymmetry state combinations
#####################################

# Map clonal states onto subclonal coordinates

obs_clonal_sorted <- obs_clonal[order(obs_clonal$coord),]
obs_subclonal_sorted <- obs_subclonal[order(obs_subclonal$coord),]
coord_clonal <- obs_clonal$coord
coord_subclonal <- obs_subclonal$coord

indices <- findInterval(
  obs_clonal_sorted$coord,
  obs_subclonal_sorted$coord
) + 1
indices[indices > length(obs_subclonal$states_subclonal)] <- length(obs_subclonal$states_subclonal)
print(head(indices,20))

AB <- data.frame(
  subclonal = obs_subclonal_sorted$states_subclonal[indices],
  clonal = obs_clonal_sorted$states_clonal
)
print(head(AB))
print(table(AB[,1]))
print(table(AB[,2]))


combined_states <- table(paste(AB$subclonal, AB$clonal, sep = ","))
Mat <- rbind(
  names(combined_states),
  as.numeric(combined_states)
)
write.table(
  Mat,
  file = file.path(output_dir, paste0(input_fn, ".CloneAS")),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)