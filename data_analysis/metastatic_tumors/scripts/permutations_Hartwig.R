library(dplyr)
library(MASS)
library(ggplot2)
library(ggdist)
library(svglite)
library(ggridges)
source("permutations_for_chroms_with_multiallelic_vars.R")

colors_list=c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")

dataAll <- read.table("../results/ALL_Hartwig_bi_multi.txt", sep="\t", header=TRUE)
dataAll$OE <- dataAll$n_multi/(((dataAll$n_bi+dataAll$n_multi)/2500000000)*(dataAll$n_bi+dataAll$n_multi))
names(dataAL)
####
dataAll<-dataAll[!duplicated(dataAll$sample),]

chrom<-read.table("../results/ALL_Hartwig_bi_multi.by_chrom.txt",sep='\t',header = T) #### your table with number of chroms
chrom$sample <- substr(chrom$sample, 1, nchar(chrom$sample) - 1)

cbind(chrom,rowSums(chrom[,c(7:14)]),rowSums(chrom[,c(3:6)]))->chrom
names(chrom)[c(15,16)]<-c('MultA', 'BiA') ### add column for number of biallelic and multiallelic mutations

merge(dataAll[,c(1,2,9,10)],chrom,by=("sample"))->CHR_OE ## keep few columns and merge with chrom data

# draw plots for all tumors with enrichment more than 10 and more than 5 multiallelic mutations
#IDs<-unique(CHR_OE$sample[CHR_OE$Multiallelic_mutations>=5 & CHR_OE$OE>=10])
IDs=c("CPCT02020393","CPCT02060274","CPCT02070070")
for (id in IDs) {
  p1<-permute_mult_sample(CHR_OE,1000,id)
  filename <- paste0("../plots/LAD/", id, "_plot.png")
  ggsave(filename=filename, plot = p1,device = "png", path = NULL,
       # scale = 1, width = 480, height = 280, units = c("mm"),
       scale = 1, width = 9, height = 10, units = c("cm"),
       dpi = 300, limitsize = TRUE)
}
