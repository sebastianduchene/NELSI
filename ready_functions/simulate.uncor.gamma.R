simulate.uncor.gamma <- function(tree, params = list(shape = 98, rate = 4361)){
    shape.gamma <- params$shape
    rate.gamma <- params$rate
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rgamma(n = length(tree$edge.length), shape = shape.gamma, rate = rate.gamma)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
}
