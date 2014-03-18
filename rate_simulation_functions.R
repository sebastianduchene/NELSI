#################
# RATE SIMULATION FUNCTIONS
#################



###############
# Constant (clock-like evolution)
#
#
###########
simulate.clock <- function(tree, rate = 0.02, noise = 0.0001){
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

#######################################################################


###############
# tdepup 1.1
# TIME DEPENDENT RATE EVOLUTION EXP DECAY
#
###########
simulate.exp.rate <- function(tree, param1 = 5, param2 = 2, noise = 0.01){
    require(phangorn)
    require(geiger)
    fun.rate <- function(x, param1 = 5, param2 = 2){
        if(any(x >= 0)){
            return((1 / (x + param1)) /param2)
        }else{
            stop("x is cannot be a negative number")
        }
    }
    #
    data.matrix <- get.tree.data.matrix(tree)
    b.times <- data.matrix[, 4]
    data.matrix[, 5] <- abs(fun.rate(b.times, param1 = param1, param2 = param2) + rnorm(nrow(data.matrix), mean = 0, sd = noise)) + 0.001
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    #
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <-c("phylogram", "tree.data.matrix")
    return(res)
}

#######################################################################

###############
# tdepup 1.2
# TIME DEPENDENT RATE EVOLUTION EXP DECAY - INTEGRATED
#
###########
simulate.exp.rate.int <- function(tree, param1 = 8, param2 = 7, noise = 0.001){
    require(phangorn)
    require(geiger)
    fun.rate <- function(x, param1 = 8, param2 = 7){
        if(any(x >= 0)){
            return((1 / (x + param1)) /param2)
        }else{
            stop("x is cannot be a negative number")
        }
    }
    # Reminder: tree.data.matrix column 5 is the rate, column 6 is the branch length substitutions, and 7 is the branch length in time.
    
    data.matrix <- get.tree.data.matrix(tree)
    b.times <- c(rep(0, length(tree$tip.label)), branching.times(tree))
    names(b.times) <- 1:length(b.times)
    # Now we take the ages of the parent nodes and the nodes.
    
    ratetemp <- vector()
    for(i in 1:length(tree$edge.length)){
    	parentage <- b.times[as.character(data.matrix[i,2])]
    	daughterage <- b.times[as.character(data.matrix[i,3])]
    	ratetemp[i] <- integrate(fun.rate, lower = daughterage, upper = parentage)$value / data.matrix[i,7]
    }
    
    data.matrix[, 5] <- abs(ratetemp + rnorm(nrow(data.matrix), mean = 0, sd = noise)) + 0.001
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    #
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <-c("phylogram", "tree.data.matrix")
    return(res)
}

#######################################################################



###############
# autocor
# AUTOCORRELATED RATE EVOLUTION
#USING THORNE ET AL. 1998 MODEL WITH HYPERPARAMETER TO CONTROL AUTOCORRELATION. NOTE THAT LOWER VALUES FOR v IMPLY HIGHER AUTOCORRELATION
###########
simulate.autocor.rate <- function(tree, initial.rate = 0.01, v = 0.3){
    require(phangorn)
    require(geiger)
    data.matrix <- get.tree.data.matrix(tree)
    while(any(is.na(data.matrix[, 5])) | any(is.nan(data.matrix[, 5]))){
        data.matrix[1, 5] <- initial.rate
        for(i in 2:nrow(data.matrix)){
            parent.node <- data.matrix[i, 2]
            preceeding.parent <- data.matrix[, 2][data.matrix[ ,3] == parent.node]
            preceeding.parent.brage <- data.matrix[, 4][data.matrix[, 2] == preceeding.parent][1]
            preceeding.parent.brrate <- data.matrix[, 5][data.matrix[, 2] == preceeding.parent][1]
            if(!(is.na(preceeding.parent.brrate)) & !(is.nan(preceeding.parent.brrate)) & (parent.node %in% data.matrix[, 3])){
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), sd = v * abs(data.matrix[i, 4] - preceeding.parent.brage)^0.5))
            }else if(!(parent.node %in% data.matrix[, 3])){
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(initial.rate)), sd = sqrt(initial.rate)))
            }
        }
    }
    data.matrix[, 6]  <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
}
#######################################################################

#############
# tdepup 1
# TIME DEPENDENT RATE EVOLUTION LOG DECAY
###########
simulate.timedep.rate.log <- function(tree, stdev = 0.001, min.rate = 0.015){
    data.matrix <- get.tree.data.matrix(tree)
    t.max <- max(data.matrix[, 4])
    b.times <- data.matrix[, 4]
#    data.matrix[, 5]  <- min.rate + ((log((abs(b.times - t.max) + 1) / t.max) + 4 ) / 200) + rnorm(length(b.times), mean = 0, sd = stdev)
	data.matrix[, 5] <- min.rate + ((-log((abs(b.times) + 1) / t.max) + 4 ) / 200) + rnorm(length(b.times), mean = 0, sd = stdev)
    data.matrix[, 6] <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
}
#######################################################################

###########
# tdepgaus
# GAUSSIAN RATE EVOLUTION
###########
simulate.gaussian.rate <- function(tree, spread = 0.3, noise.sd = 0.0008){
    data.matrix <- get.tree.data.matrix(tree)
    mean.gauss <- max(data.matrix[, 4]) / 2 # use to set the peak of the function at the median median branch time
    #mean.gauss <- max(allnode.times(tree)) / 2 # Use to set the peak of the function at the mean height of the tree
    spread.gauss <- spread * mean.gauss
    data.matrix[, 5] <- 10*abs(0.1 * dnorm(data.matrix[, 4], mean = mean.gauss, sd = spread.gauss) + rnorm(nrow(data.matrix), mean = 0, sd = noise.sd) + 0.002)
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
}
##############################################################

############
# uncorlog
# UNCORRELATED LOGNORMAL
###########
simulate.lognormal.uncor.rate <- function(tree, mean.log = -3.9, sd.log = 0.1){
	data.matrix <- get.tree.data.matrix(tree)
	branch.rates <- rlnorm(n = length(tree$edge.length), meanlog = mean.log, sdlog = sd.log)
	data.matrix[, 5] <- branch.rates
	data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
	tree$edge.length <- data.matrix[, 6]
	res <- list(tree, data.matrix)
	names(res) <- c("phylogram", "tree.data.matrix")
	return(res)
}
##############################################################

############
# tdepdown
# TIME DEPENDENT SIGMOIDAL
###########
simulate.timedep.rate.logit <- function(tree, stdev = 0.005, min.rate = 0.025){
    data.matrix <- get.tree.data.matrix(tree)
    t.max <- max(data.matrix[, 4])
    b.times <- data.matrix[, 4]
    data.matrix[, 5] <- abs(min.rate + (1 / (1 + ((1.15 ^ t.max) * exp(-(b.times/3.2))))) / ((t.max/2.8)*3) + rnorm(n = length(b.times), mean = 0, sd = stdev))
    data.matrix[, 6] <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
}
##############################################################
