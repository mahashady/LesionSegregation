simulate_controls_duration<-function (tree, node, interdist, iterations) {
  simul_dist <- NULL
  iter = 0
  while (iter < iterations){
    random_node = sample((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)), 1)
    node_ancestors <- getAncestors(tree, random_node)
    node_children = get_node_children(random_node, tree)
    node_all_children = get_all_node_children(random_node, tree)
    if ((length(node_ancestors) >= interdist + 1) & length(node_children) > 2 & random_node != node & nodeheight(tree, node) > 60){
      print(random_node)
      ancestor_node = node_ancestors[interdist]  
      dist <- nodeheight(tree, random_node) - nodeheight(tree, ancestor_node)
      simul_dist<-c(simul_dist, dist)
      iter <- iter + 1
    }
  }
  return(simul_dist)
}

simulate_controls_leaves<-function (tree, node, allowed_distance, iterations) {
  simul_number_children <- NULL
  node_time = nodeheight(tree, node)
  close_nodes_children = NULL
  for (random_node in (Ntip(tree)+2):(Ntip(tree)+Nnode(tree))){
    random_node_time = nodeheight(tree, random_node)
    random_node_ancestors <- getAncestors(tree, random_node)
    random_node_children = get_node_children(random_node, tree)
    if (abs(random_node_time-node_time) <= allowed_distance & random_node != node & length(random_node_ancestors)>=2 & length(random_node_children) >=2 & random_node_time > 60){
      n_children = length(get_all_node_leaves(random_node, tree))
      close_nodes_children <- c(close_nodes_children, n_children)
    }
  }
  print("close nodes_children")
  print(close_nodes_children)
  if (length(close_nodes_children) != 0){
    simul_number_children <- sample(close_nodes_children, iterations, replace=TRUE)
  }else{
    simul_number_children <- NULL
  }
  print(simul_number_children)
  return(simul_number_children)
}


simulate_controls_branching_end_node<-function (tree, node, time_after, iterations) {
  simul_branchings <- NULL
  iter = 0
  while (iter < iterations){
    random_node = sample((Ntip(tree)+1):c(Ntip(tree)+Nnode(tree)), 1)
    node_ancestors <- getAncestors(tree, random_node)
    node_children = get_node_children(random_node, tree)
    check_longest_offspring = FALSE
    all_random_node_leaves = get_all_node_leaves(random_node, tree)
    for (leaf in all_random_node_leaves){
      if (nodeheight(tree,leaf)-nodeheight(tree,random_node)>time_after){
        check_longest_offspring = TRUE
      }
    }
    if ((length(node_ancestors) >= 2) & length(node_children) >= 2 & random_node != node & check_longest_offspring & nodeheight(tree, random_node)>60){
      n_branchings<-length(get_branching_points_in_time(tree, random_node, time_after))
      print(n_branchings)
      simul_branchings<-c(simul_branchings, n_branchings)
      iter <- iter + 1
    }
  }
  return(simul_branchings)  
}

simulate_controls_branching_start_node<-function (tree, node, time_after, iterations) {
  simul_branchings <- NULL
  iter = 0
  while (iter < iterations){
    random_node = sample((Ntip(tree)+1):c(Ntip(tree)+Nnode(tree)), 1)
    node_ancestors <- getAncestors(tree, random_node)
    node_children = get_node_children(random_node, tree)
    check_longest_offspring = FALSE
    all_random_node_leaves = get_all_node_leaves(random_node, tree)
    for (leaf in all_random_node_leaves){
      if (nodeheight(tree,leaf)-nodeheight(tree,random_node)>time_after){
        check_longest_offspring = TRUE
      }
    }
    if ((length(node_ancestors) >= 1) & length(all_random_node_leaves) >= 3 & random_node != node & check_longest_offspring & nodeheight(tree, random_node)>60){
      n_branchings<-length(get_branching_points_in_time(tree, random_node, time_after))
      print(n_branchings)
      simul_branchings<-c(simul_branchings, n_branchings)
      iter <- iter + 1
    }
  }
  return(simul_branchings)  
}

simulate_controls_drivers<-function (tree, node, drivers, iterations) {
  simul_drivers <- NULL
  iter = 0
  while (iter < iterations){
    random_node = sample((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)), 1)
    node_ancestors <- getAncestors(tree, random_node)
    node_children = get_node_children(random_node, tree)
    if (length(node_ancestors) >= 1 & length(node_children) >= 2 & random_node != node & nodeheight(tree, node) > 60){
      print(random_node)
      simul_drivers<-c(simul_drivers, sum(random_node %in% drivers))
      iter <- iter + 1
    }
  }
  return(simul_drivers)
}

simulate_controls_drivers_up<-function (tree, node, drivers, dist2root, iterations) {
  simul_drivers <- NULL
  iter = 0
  while (iter < iterations){
    random_node = sample((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)), 1)
    node_ancestors <- getAncestors(tree, random_node)
    node_children = get_node_children(random_node, tree)
    if (length(node_ancestors) >= 1 & length(node_children) >= 2 & random_node != node & nodeheight(tree, node) > 60 & length(node_ancestors) >= dist2root){
      print(random_node)
      drivers_up = ifelse(length(intersect(node_ancestors[1:dist2root], drivers)) >=1, 1,0)
      simul_drivers<-c(simul_drivers, drivers_up)
      iter <- iter + 1
    }
  }
  return(simul_drivers)
}
