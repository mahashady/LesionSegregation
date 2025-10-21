#### Function to number of homologous chromosomes for each chromosome after N divisions. It starts with vector that have 2, 22 times
simulate_divisions <- function(input_vector, N) {
  if (N == 0) {return(input_vector)}
  for (iteration in 1:N) {
    for (i in 1:length(input_vector)) {
      current_value <- input_vector[i]
      
      if (current_value == 2) {
        probabilities <- c(1/4, 1/2, 1/4)
        modified_value <- sample(c(2, 1, 0), size = 1, prob = probabilities)
      } else if (current_value == 1) {
        probabilities <- c(1/2, 1/2)
        modified_value <- sample(c(1, 0), size = 1, prob = probabilities)
      } else {
        modified_value <- 0
      }
      
      input_vector[i] <- modified_value
    }
  }
  if (sum(input_vector)==0){input_vector[1]=1}
  return(input_vector)
}


permute_mult_sample <- function(data, iter, id) {
  iv<-rep(2,22)
  data[data$sample==id, ]->local_data
  print(head(local_data))
  dist=data.frame();
  
  for (mrca in 0:4){
    mrca1=mrca+1;
    for (i in 1:iter){
      #  data=CHR_OE_M;column1="MultA";group_column="sample";group='CPCT02010122';
      
      Obs=0;EXP=0
      
      #group="CPCT02410017"
      result_vector <- simulate_divisions(iv, mrca);
      local_data$BiA/local_data$Total_mutations*(result_vector/sum(result_vector))->probabilities; #these are probabilities of multiallelic sites on each chromosome depending on bi-allelic mutation rate and number of homologous chroms with multi sites
      # Select a box based on the probabilities - selected chromosomes with multi-sites based on the probabilities above 
      selected_box <- sample(1:22, size = sum(local_data$MultA), prob = probabilities, replace = TRUE)
      #print (c(group,sum(CHR_OE_M$MultA[CHR_OE_M$sample==group]),selected_box,
      #       sum(as.numeric(CHR_OE_M$MultA[CHR_OE_M$sample==group]>0.1)),
      #       length(unique(selected_box))));
      EXP<-length(unique(selected_box))
      Obs<-sum(as.numeric(local_data$MultA>0.1))
      rbind(dist,c(EXP,Obs,mrca1))->dist}
    names(dist)<-c("exp","obs","mrca");}
  #p1 <-ggplot(dist, aes(exp)) +
  # geom_histogram(binwidth = 1, color = "black", fill = "white") +
  #  facet_wrap(~mrca, strip.position = "top", scales = "fixed",ncol = 5) +
  #  xlab("#chr with multiallelic mutations") +
  #  ylab("density")+
  #  theme_classic() +
  # Add red arrow for each panel
  #  geom_segment(
  #    aes(x = obs, y = Inf, xend = obs, yend = 0),  # Set yend to Inf
  #    color = "red", arrow = arrow(),  size = 1,
  #    data = data.frame(mrca = unique(dist$mrca), obs = unique(dist$obs)),
  #    inherit.aes = FALSE) +
  #  stat_density(aes(y = after_stat(density)), geom = "line", position = "identity", color = "black")+
  #  ggtitle(paste(id)) +
  #  theme(plot.title = element_text(hjust = 0.5))  # Center the title horizontally)
  n_bins=max(dist$exp)
  print(n_bins)
  p1 <- (ggplot(dist, aes(x=exp,y=factor(mrca, levels=c(5,4,3,2,1)))) 
         + geom_density_ridges(alpha=0.6,scale=0.8,stat = "binline", bins=n_bins,fill="grey", col="darkgrey")
         + geom_vline(xintercept = Obs, size=1)
#         + annotate(geom="text", label="observed", x=Obs-5, y=0.6, size=5)
         + theme(plot.title = element_text(hjust = 0.5))
         + theme_bw()
         + scale_fill_manual(values=colors_list[2:6])
         + scale_color_manual(values=colors_list[2:6])
         + ylab("LAD")
         + theme(legend.position = "none")
         + xlab("N chroms with multiallelic sites")
         + ggtitle(id)
  )
  return(p1)
  
}
