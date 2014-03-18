
simulate.tdep.ho <- function(tree, params = list(mu = 0.015, srate = 0.035, lambda = 0.1, noise = 0.001)){
    require(phangorn)
    require(geiger)
    mu <- params$mu
    srate <- params$srate
    lambda <- params$lambda
    noise <- params$noise
    fun.rate <- function(x, m = mu, s = srate, lam = lambda){
        if(any(x >= 0)){
            return(s + (m * exp(-lam * x)))
        }else{
            stop("x is cannot be a negative number")
        }
    }

    data.matrix <- get.tree.data.matrix(tree)
    b.times <- c(rep(0, length(tree$tip.label)), branching.times(tree))
    names(b.times) <- 1:length(b.times)

    ratetemp <- vector()
    for(i in 1:length(tree$edge.length)){
    	parentage <- b.times[as.character(data.matrix[i,2])]
    	daughterage <- b.times[as.character(data.matrix[i,3])]
    	ratetemp[i] <- integrate(fun.rate, lower = daughterage, upper = parentage)$value / data.matrix[i,7]
    }

    data.matrix[, 5] <- abs(ratetemp + rnorm(nrow(data.matrix), mean = 0, sd = noise))
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]

    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <-c("phylogram", "tree.data.matrix")
    return(res)
}
