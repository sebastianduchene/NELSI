simulate.white.noise <-
function(tree, params = NULL){
    data.matrix <- get.tree.data.matrix(tree)
    means.rates <- data.matrix[, 7] / sum(data.matrix[, 7])
    branch.rates <- sapply(means.rates, function(x) log(rlnorm(1, x, x)))
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, tree.data.matrix = data.matrix)
    class(res) <- "ratesim"
    return(res)
}
