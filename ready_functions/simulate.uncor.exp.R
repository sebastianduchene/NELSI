simulate.uncor.exp <- function(tree, params = list(mean.exp = 0.001)){
    mean.exp <- params$mean.exp
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rexp(n = length(tree$edge.length), rate = 1 / mean.exp)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
