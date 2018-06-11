get.external.branch.length <- function(tr){
    tr$edge.length[!tr$edge[, 2] %in% tr$edge[, 1]]
}
