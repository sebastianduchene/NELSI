simulate.uncor.gamma <-
function(tree, params = list(mean = 0.007, shape = 3.98, rate = NULL)){
    if(is.null(params$rate)) params$rate <- params$shape/params$mean
    rate.gamma <- params$rate
    shape.gamma <- params$shape
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rgamma(n = length(tree$edge.length), shape = shape.gamma, rate = rate.gamma)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
