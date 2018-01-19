get.internal.branch.length <- function(tr){
    internal_branches <- which(tr$edge[, 2] %in% tr$edge[, 1])
    return(tr$edge.length[internal_branches])
}
