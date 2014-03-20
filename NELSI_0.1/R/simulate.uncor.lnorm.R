simulate.uncor.lnorm <-
function(tree, params = list(mean.log = -3.9, sd.log = 0.1)){
    mean.log <- params$mean.log
    sd.log <- params$sd.log
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rlnorm(n = length(tree$edge.length), meanlog = mean.log, sdlog = sd.log)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
