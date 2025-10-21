args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]

fnf<-paste("/workspace/projects/lesion_segregation/mice/results/HMM/input_HMM/", fn,sep = "")

library(HMM)
library(betareg)
library(flexmix)
library(mixtools)

data<-read.table(file=fnf, header=F, sep=" ")####
colnames(data) <- c("Chr", "Pos", "Mutation", "Context", "Cov_ref", "Cov_alt", "Sample_id", "bi_multi", "Substitution_type")

data$Substitution_type<-as.character(data$Substitution_type)
data$bi_multi<-as.character(data$bi_multi)
data$Chr<-as.character(data$Chr)

SM <- data[data$bi_multi == 'B' & (data$Substitution_type == 'A_N' | data$Substitution_type == 'T_N') & data$Chr != 'X',]
SM$Chr<-as.numeric(SM$Chr)
SM$Cov_alt<-as.numeric(SM$Cov_alt)
SM$vaf <- SM$Cov_alt/(SM$Cov_alt+SM$Cov_ref)
SM$coord <- as.numeric(SM$Chr)*1000000000+SM$Pos

fit <- normalmixEM(SM$vaf, k = 2) #fit vaf distribution as a mixture of two normal distributions
sum(fit$posterior[,2]>0.8)
sum(fit$posterior[,1]>0.8)


TM<-matrix(c(rep(0.002,9)),3) 
diag(TM) <- rep(0.192,3)

#run HMMas on the whole set of mutations
hmm=initHMM(c("A0","A1","A2" ),c("T_N","A_N"),
            transProbs=TM, 
            emissionProbs=rbind(c(.9,.1),
                                c(.5,.5),
                                c(.1,.9))) 

#run Baum-Welch and Viterbi for clonal mutations
bw2=baumWelch(hmm,SM[(SM$Substitution_type == 'A_N' | SM$Substitution_type == 'T_N') & fit$posterior[,2] > 0.8,]$Substitution_type, 15)
viterbi(bw2$hmm,SM[(SM$Substitution_type == 'A_N' | SM$Substitution_type == 'T_N') & fit$posterior[,2] > 0.8,]$Substitution_type) -> l2

#run Baum-Welch and Viterbi for subclonal mutations
bw1=baumWelch(hmm,SM[(SM$Substitution_type == 'A_N' | SM$Substitution_type == 'T_N') & fit$posterior[,1] > 0.8,]$Substitution_type, 15)
viterbi(bw1$hmm,SM[(SM$Substitution_type == 'A_N' | SM$Substitution_type == 'T_N') & fit$posterior[,1] > 0.8,]$Substitution_type) -> l1

#create matrix with emissions and number of sites in each state separately for clonal and subclonal mutations and write to the output
EMP <- as.vector(c('Emissions_bw1_bw2',bw1$hmm$emissionProbs[,1],bw2$hmm$emissionProbs[,1]))
STP <- c('Statecounts_bw1_bw2',sum(as.numeric(l1=='A0')),sum(as.numeric(l1=='A1')),sum(as.numeric(l1=='A2')),
       sum(as.numeric(l2=='A0')),sum(as.numeric(l2=='A1')),sum(as.numeric(l2=='A2')))
Bw_matrix <- rbind(EMP,STP)
fnf1 <- paste("/workspace/projects/lesion_segregation/mice/results/HMM/HMM_ploidy/", fn,".BW",sep = "")
write.table((Bw_matrix), file=fnf1, sep=" ", row.names=FALSE, quote = FALSE, col.names = FALSE)

#output boundary vafs for clonal and subclonal mutations, median vaf of clonal and subclonal mutations, number of clonal and subclonal mutations
MAX <- min(max(SM[fit$posterior[,1]>0.8,]$vaf), max(SM[fit$posterior[,2]>0.8,]$vaf))
MIN<-max(min(SM[fit$posterior[,1]>0.8,10]),min(SM[fit$posterior[,2]>0.8,10]))
Distr_properies<-c('vafMed_vafboundaries_Distr_size',median(SM[fit$posterior[,1]>0.5,10]),median(SM[fit$posterior[,2]>0.5,10]),
MAX,MIN,sum(as.numeric(fit$posterior[,1]>0.8)),sum(as.numeric(fit$posterior[,2]>0.8)))
fnf2<-paste("/workspace/projects/lesion_segregation/mice/results/HMM/HMM_ploidy/", fn,".Clonesize",sep = "")
write.table((Distr_properies),file=fnf2,sep=" ",row.names=FALSE,quote = FALSE,col.names = FALSE)

#output the asymmetry states of clonal and subclonal mutations
indices <- findInterval(SM[fit$posterior[,2]>0.8,]$coord, sort(SM[fit$posterior[,1]>0.8,]$coord))+1
indices[indices > length(l1)] <- length(l1)
AB <- as.data.frame(cbind(l1[indices],l2))
combined_states <- table(paste(AB[,1], AB[,2], sep = ","))
as.numeric(t(as.matrix(combined_states)))->CS
as.vector(names(combined_states))->nam
rbind(nam,as.numeric(CS))->Mat
fnf3<-paste("/workspace/projects/lesion_segregation/mice/results/HMM/HMM_ploidy/", fn,".CloneAS",sep = "")
write.table((Mat),file=fnf3,sep=" ",row.names=FALSE,quote = FALSE,col.names = FALSE)