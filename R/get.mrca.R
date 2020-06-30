get.mrca <- function(tr, tips){
  root <- tr$edge[, 1][!tr$edge[, 1] %in% tr$edge[, 2]][1]
  paths_to_root <- list()
  for(i in 1:length(tips)){
    nodes_in_path <- tips[i]
    ancestor <- 0
    while(ancestor != root){
      ancestor <- tr$edge[nodes_in_path[length(nodes_in_path)] == tr$edge[, 2], 1]
      nodes_in_path <- c(nodes_in_path, ancestor)
    }
    paths_to_root[[i]] <- nodes_in_path
  }
  rintersect <- function(x){
    if(length(x) == 2){
      return(intersect(x[[1]], x[[2]]))
      }else{
        return(intersect(x[[1]], rintersect(x[-1])))
      }
  }
  if(length(tips) > 1){
    mrca <- rintersect(paths_to_root)[1]
    }else{
      mrca <- paths_to_root[[1]][2]
    }
  return(mrca)
}


