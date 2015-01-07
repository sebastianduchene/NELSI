simulate.dpp <-
function(tree, params = list(alpha = 1, shape = 10, rate = 0.0001)){
    shape.gamma <- params$shape
    rate.gamma <- params$rate
    alpha <- params$alpha
    data.matrix <- get.tree.data.matrix(tree)
    nbranches <- nrow(tree$edge)
    cats <- sample(1:nbranches, 1, prob = rdirichlet(1, rep(alpha, nbranches)))
    branch.cats <- sample(1:cats, nbranches, replace = T, prob = rdirichlet(1, rep(alpha, cats)))
    branch.rates <- rgamma(n = cats, shape = shape.gamma, rate = rate.gamma)
    names(branch.rates) <- 1:cats
    branch.rates <- branch.rates[branch.cats]
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
