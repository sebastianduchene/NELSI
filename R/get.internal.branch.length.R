get.internal.branch.length <- function(tr){
    return(which(tr$edge[, 2] %in% tr$edge[, 1]))
}
