cran_packages=c("stringr","ape","remotes")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}
options(stringsAsFactors = FALSE)


muts_file_initial <- "/workspace/projects/lesion_segregation/blood/data/annotated_mut_set_PX001_2_01_standard_rho01"
muts_file_reassigned <- "/workspace/projects/lesion_segregation/blood/data/annotated_mut_set_PX001_2_01_standard_rho01.reassigned"

load(muts_file_initial)
muts_initial=filtered_muts$COMB_mats.tree.build$mat

load(muts_file_reassigned)
muts_reassigned=filtered_muts$COMB_mats.tree.build$mat

print(nrow(muts_initial))
print(nrow(muts_reassigned))

result <- merge(muts_initial[,c("mut_ref", "node")], muts_reassigned[,c("mut_ref", "node")], by= "mut_ref", suffixes=c(".initial", ".reassigned"))
print(head(result))
print(nrow(result[result$node.initial != result$node.reassigned]))
result[result$mut_ref == "17-58740698-A-T",]