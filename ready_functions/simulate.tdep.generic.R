simulate.tdep.generic <- function(tree, noise = 0.001, fun = function(time){ 0.035 * (0.015 * exp(-0.1 * time)) }){
    require(phangorn)
    require(geiger)
    fun.rate <- function(time){
        if(any(time >= 0)){
            return(fun(time))
        }else{
            stop("x is cannot be a negative number")
        }
    }

    data.matrix <- get.tree.data.matrix(tree)
    node.ages <- allnode.times(tree)
    b.times <- c(rep(0, length(tree$tip.label)), node.ages[(length(tree$tip.labels)+1):length(node.ages)])
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
