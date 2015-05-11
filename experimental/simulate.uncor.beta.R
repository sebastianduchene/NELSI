simulate.uncor.beta <-
function(tree, params = list(shape1 = 0.4, shape2 = 0.4, scale_param = 150, center_param = 0.0015)){
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- (rbeta(n = length(tree$edge.length), shape1 = params$shape1, shape2 = params$shape2) / params$scale_param) + params$center_param
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
