library(ggplot2)
library(HMM)
library(this.path)

args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]

setwd(dirname(this.path()))
print(getwd())
fnf <- paste("../results/HMM/input_HMM/", fn,sep = "")

data <- read.table(file=fnf, header=F, sep=" ")
colnames(data) <- c("Chr", "Pos", "Mutation", "Context", "Cov_ref", "Cov_alt", "Sample_id", "bi_multi", "Substitution_type")
print(head(data))

TM <- matrix(c(rep(0.002,25)),5) 
diag(TM) <- rep(0.192,5)


hmm=initHMM(c("A0","A1","A2","A3","A4" ),c('T_N','A_N'),
            transProbs=TM, 
            emissionProbs=rbind(c(.9,.1),
                                c(.7,.3),
                                c(.5,.5),
                                c(.3,.7),
                                c(.1,.9))) 




data$Substitution_type <- as.character(data$Substitution_type)
data$bi_multi <- as.character(data$bi_multi)
data$Chr <- as.character(data$Chr)
data[data$bi_multi=='B',]$Cov_alt <- as.integer(as.character(data[data$bi_multi=='B',]$Cov_alt))


#RUN HMMas (HMM for WC asymmetry)
print("Running HMMas")
bw <- baumWelch(hmm,data[data$Substitution_type == 'A_N' | data$Substitution_type == 'T_N',]$Substitution_type,25,delta=1E-4, pseudoCount=0)
viterbi(bw$hmm,data[data$Substitution_type == 'A_N' | data$Substitution_type == 'T_N',]$Substitution_type) -> l
  

TM<-matrix(c(rep(0.002,9)),3) 
diag(TM) <- rep(0.3293333,3)


hmm=initHMM(c("A0","A1","A2" ),c("M","B"),
            transProbs=TM, 
            emissionProbs=rbind(c(.02,.98),
                                c(.08,.92),
                                c(.16,.84))) 

#RUN HMMmulti 
print("Running HMMmulti")
bwM <- baumWelch(hmm,data$bi_multi,25)
viterbi(bwM$hmm,data$bi_multi) -> lm


#####
SM <- data[data$Substitution_type == 'A_N' | data$Substitution_type == 'T_N',]
SM$Cov_alt <- as.integer(as.character(SM$Cov_alt))
SM$vaf <- SM$Cov_alt/(SM$Cov_alt+SM$Cov_ref) 
st='A0'; tr=.15

#Calculate WC asymmetry for subclonal mutations in HMMas == "A0" (vaf cut-off < 0.15)
SubAS <- sum(as.numeric(SM[l == st & SM$bi_multi == 'B' & SM$Chr != 'X' & SM$vaf < tr,]$Substitution_type == 'A_N'))/
sum(as.numeric(SM[l == st & SM$bi_multi == 'B' & SM$Chr != 'X' & SM$vaf < tr,]$Substitution_type == 'T_N'))

#Calculate WC asymmetry for clonal mutations in HMMas == "A0" (vaf cut-off > 0.25)
tr=.25
ClonAS<-sum(as.numeric(SM[l == st & SM$bi_multi == 'B' & SM$Chr != 'X' & SM$vaf > tr,]$Substitution_type == 'A_N'))/
sum(as.numeric(SM[l == st & SM$bi_multi == 'B' & SM$Chr != 'X' & SM$vaf > tr,]$Substitution_type == 'T_N'))


data[data$bi_multi == 'B' & data$Chr != 'X',] -> Aut
data[data$bi_multi == 'B' & data$Chr == 'X',] -> SEX

as.integer(as.character(SEX$Cov_alt)) -> SEX$Cov_alt
as.integer(as.character(Aut$Cov_alt)) -> Aut$Cov_alt

SE <- median(SEX$Cov_alt/(SEX$Cov_alt + SEX$Cov_ref))
Au <- median(Aut$Cov_alt/(Aut$Cov_alt + Aut$Cov_alt))

out <- c(fn, ClonAS, SubAS, Au, SE, sum(as.numeric(l=='A0')), sum(as.numeric(l=='A1')), sum(as.numeric(l=='A2')),
sum(as.numeric(l=='A3')), sum(as.numeric(l=='A4')),
sum(as.numeric(lm=='A0')), sum(as.numeric(lm=='A1')), sum(as.numeric(lm=='A2')))


fnf1 <- paste("../results/HMM/output_HMM/", fn,".summary",sep = "")
fnf2 <- paste("../results/HMM/output_HMM/", fn,".ASHMM",sep = "")
fnf3 <- paste("../results/HMM/output_HMM/", fn,".PloidyHMM",sep = "")
fnf4 <- paste("../results/HMM/output_HMM/", fn,".bwASEM",sep = "") ### bw$hmm$emissionProbs
fnf5 <- paste("../results/HMM/output_HMM/", fn,".bwPlEM",sep = "") ### bwM$hmm$emissionProbs


write.table(t(out), file = fnf1, sep=" ", row.names=FALSE, quote = FALSE, col.names = FALSE) #output sample summary
write.table((l), file = fnf2, sep=" ", row.names=FALSE, quote = FALSE, col.names = FALSE) #output HMMas state per site
write.table((lm),file = fnf3,sep=" ",row.names=FALSE,quote = FALSE,col.names = FALSE) #output HMMmulti state per site
write.table((bw$hmm$emissionProbs), file=fnf4, sep=" ", row.names=FALSE, quote = FALSE, col.names = FALSE) #output emission probabilities of HMMas
write.table((bwM$hmm$emissionProbs),file=fnf5,sep=" ",row.names=FALSE,quote = FALSE,col.names = FALSE) #output emission probabilities of HMMmulti