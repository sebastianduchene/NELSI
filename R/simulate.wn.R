simulate.wn <- 
function(tree, params = list(rate = 0.006, var = 0.000001)){
    data.matrix <- get.tree.data.matrix(tree)
    clocksubst <- tree$edge.length * params[[1]]
    clocksubstscaled <- as.numeric(scale(clocksubst)) + 1
    for(i in 1:length(clocksubst)) data.matrix[, 6][i] <- rlnorm(1, log(clocksubst[i]), abs(params[[2]] * clocksubstscaled[i]))
    branch.rates <- data.matrix[, 6] / data.matrix[, 7]
    data.matrix[, 5] <- branch.rates
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, tree.data.matrix = data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}