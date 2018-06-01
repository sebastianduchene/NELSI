get.deepest.br.length <- function(tr){
  root <- tr$edge[!tr$edge[, 1] %in% tr$edge[, 2], 1][1]
  descending.branches <- tr$edge.length[tr$edge[, 1] == root]
  sum(descending.branches)
}
