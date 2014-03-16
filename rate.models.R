# This script contains the functions to simulate under different rate models

require(ape)
require(geiger)
require(phangorn)
source("utilities.R")

######
#Testing

tr1 <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 20)
s1 <- simulate.clock(tr1)

params <- list(rate = 0.002, noise = 0.00001)

s2 <- simulate.clock(tr1, params)
######



# Evolution under a strict clock
simulate.clock <- function(tree, params = list(rate = 0.02, noise = 0.0001)){
    rate <- params$rate
    noise <- params$noise
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rep(rate, times = length(tree$edge.length))
    branch.rates <- abs(branch.rates + rnorm(length(tree$edge.length), mean = 0, sd = noise))
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
}
