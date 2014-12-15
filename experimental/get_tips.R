get.tips <- function(tr, node){
   child <- tr$edge[, 2][tr$edge[, 1] == node]
   if(all(child <= length(tr$tip.label))){
      return(child)
   }else{
      nodes <- child[which(child > length(tr$tip.label))]
      children <- unlist(c(child[which(child <= length(tr$tip.label))], sapply(nodes, function(x) get_tips(tr, x))))
      return(children)
   }
}
