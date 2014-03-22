get.tree.data.matrix <-
function(phylo){
    require(phangorn)
    require(geiger)
    data.matrix <- matrix(data = NA, ncol = 7, nrow = length(phylo$edge.length))
    colnames(data.matrix) <- c("branch.index", "parent.node", "daughter.node", "branch.midage", "branch.rate", "length.subst", "length.time")
    data.matrix[, 1] <- 1:length(phylo$edge.length)
    data.matrix[, 2] <- phylo$edge[ ,1]
    data.matrix[, 3] <- phylo$edge[ ,2]
    data.matrix[, 4] <- mid.edge.ages(phylo)
    data.matrix[, 7] <- phylo$edge.length
    class(data.matrix) <- "tree.data.matrix"
    return(data.matrix)
}
