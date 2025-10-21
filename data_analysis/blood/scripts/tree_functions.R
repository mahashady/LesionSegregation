get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}

get_all_node_leaves=function(node,tree){
  leaves = intersect(1:Ntip(tree), get_all_node_children(node, tree))
  return(leaves)
}

get_all_node_branchings=function(node,tree){
  branchings = intersect((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)), get_all_node_children(node, tree))
  return(branchings)
}

getAncestors=function (tree, node, type = c("all", "parent")) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  type <- type[1]
  if (type == "all") {
    aa <- vector()
    rt <- Ntip(tree) + 1
    currnode <- node
    while (currnode != rt) {
      currnode <- getAncestors(tree, currnode, "parent")
      aa <- c(aa, currnode)
    }
    return(aa)
  }
  else if (type == "parent") {
    aa <- tree$edge[which(tree$edge[, 2] == node), 1]
    return(aa)
  }
  else stop("do not recognize type")
}


nodeheight=function (tree, node, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (node == (Ntip(tree) + 1)) 
    h <- 0
  else {
    a <- setdiff(c(getAncestors(tree, node), node), Ntip(tree) + 1)
    h <- sum(tree$edge.length[sapply(a, function(x, e) which(e == x), e = tree$edge[, 2])])
  }
  h + ROOT
}


get_branching_points_in_time<-function(tree, node, time){
  branching_in_time =NULL 
  all_node_branchings = get_all_node_branchings(node, tree)
  for (branching in all_node_branchings){
    if (nodeheight(tree, branching) - nodeheight(tree, node) <= time){
      branching_in_time <- c(branching_in_time, branching)
    }
  }
  return(branching_in_time)
}
