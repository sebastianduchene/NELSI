get.external.branch.length <- function(tr){
  #set.seed(1234)
  #tr <- rtree(5)
  #plot(tr, show.tip.label = F)
  #edgelabels(round(tr$edge.length, 2)) 
  #tiplabels()
  #nodelabels()
  external.nodes <- tr$edge[!tr$edge[, 2] %in% tr$edge[, 1], 2]
  external.edges <- which(tr$edge[, 2] %in% external.nodes)
  external.edges.lengths <- tr$edge.length[external.edges]
  return(external.edges.lengths)
}