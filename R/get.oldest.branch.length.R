get.oldest.branch.length <- function(tr){
    root_node <- tr$edge[!(tr$edge[, 1] %in% tr$edge[, 2]), 1][1]
    min(tr$edge.length[tr$edge[, 1] == root_node])
}
