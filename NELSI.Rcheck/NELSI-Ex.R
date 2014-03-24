pkgname <- "NELSI"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('NELSI')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("NELSI_0.2-package")
### * NELSI_0.2-package

flush(stderr()); flush(stdout())

### Name: NELSI-package
### Title: NELSI: Nucleotide EvoLution Simulator
### Aliases: NELSI-package NELSI
### Keywords: package

### ** Examples

set.seed(1234525)

myTree <- rcoal(50)

# Simulate uncorrelated rates with default parameters:
rateTree.default <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot(rateTree.default, col.lineages = rainbow(50))



cleanEx()
nameEx("allnode.times")
### * allnode.times

flush(stderr()); flush(stdout())

### Name: allnode.times
### Title: allnode.times
### Aliases: allnode.times
### Keywords: chronogram

### ** Examples

set.seed(12345)
myTree <- rtree(10)
plot(myTree)
# See the numbering of internal nodes and tips. Note that the tip labels and the actual numbering of the tips are different.
nodelabels()
tiplabels()
allnode.times(myTree)

# Plot the tree and add the ages of the tips and internal nodes.
plot(myTree, show.tip.label = FALSE)
allTimes <- allnode.times(myTree)
allTimes
tiplabels(round(allTimes[1:10]))
nodelabels(round(allTimes[11:19]))


## The function is currently defined as
function (phylo, tipsonly = FALSE) 
{
    require(phangorn)
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 
        2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == 
        root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == 
        root.phylo, ]
    if (tipsonly == TRUE) {
        node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == 
            root.phylo, 1:length(phylo$tip.label)]
    }
    return(node.times)
  }



cleanEx()
nameEx("get.clock.data")
### * get.clock.data

flush(stderr()); flush(stdout())

### Name: get.clock.data
### Title: get.clock.data
### Aliases: get.clock.data
### Keywords: molecular-clock

### ** Examples

set.seed(12345)
myTree <- rtree(10) # Note that this tree is not ultrametric.
myTreeTimes <- allnode.times(myTree)


plot(myTree, show.tip.label = FALSE)
tiplabels(round(myTreeTimes[1:10], 2))
nodelabels(round(myTreeTimes[11:19], 2))

# Simulate rates along the tree with the uncorrlated lognormal model with default settings.
rateTree <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)

clockDataTree <- get.clock.data(rateTree, pch = 20, col = "blue")

# Linear regression for substitutions as a function of time

lmRate <- lm(substitutions ~ times, data = clockDataTree)
summary(lmRate)


## The function is currently defined as
function (rate.sim.object, tipsonly = T, ...) 
{
    phylogram <- rate.sim.object$phylogram
    chrono <- rate.sim.object$phylogram
    chrono$edge.length <- rate.sim.object[[2]][, 7]
    times <- allnode.times(chrono, tipsonly)
    substitutions <- allnode.times(phylogram, tipsonly)
    plot(times, substitutions, ...)
    return(data.frame(times, substitutions))
  }



cleanEx()
nameEx("get.lineage.time.rate")
### * get.lineage.time.rate

flush(stderr()); flush(stdout())

### Name: get.lineage.time.rate
### Title: get.lineage.time.rate
### Aliases: get.lineage.time.rate
### Keywords: substitution rate phylo

### ** Examples

set.seed(123425)
myTree <- rcoal(10)
rateTree <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot(rateTree, col.lineages = rainbow(10))

# Get the rate for the lineages ascending from the first taxon
# Find the name of the first taxon
myTree$tip.label

get.lineage.time.rate(taxon = 1, sim.rate.object = rateTree)



## The function is currently defined as
function (taxon, sim.rate.object) 
{
    tree.data.matrix <- sim.rate.object[[2]]
    chrono <- sim.rate.object[[1]]
    chrono$edge.length <- tree.data.matrix[, 7]
    taxon.init <- taxon
    if (taxon %in% tree.data.matrix[, 3]) {
        data.matrix <- tree.data.matrix
        branch.times <- vector()
        rate.time <- vector()
        repeat {
            parent <- data.matrix[, 2][data.matrix[, 3] == taxon]
            time.br <- data.matrix[, 4][data.matrix[, 3] == taxon]
            rate.br <- data.matrix[, 5][data.matrix[, 3] == taxon]
            rate.time <- c(rate.time, rate.br)
            branch.times <- c(branch.times, time.br)
            taxon <- parent
            if (!(parent %in% data.matrix[, 3])) {
                break
            }
        }
        first.rate <- rate.time[length(rate.time)]
        last.rate <- rate.time[1]
        rate.time <- c(last.rate, rate.time, first.rate)
        node.times <- allnode.times(chrono)
        root.age <- max(node.times)
        branch.times <- c(node.times[taxon.init], branch.times, 
            root.age)
        return(data.frame(branch.times, rate.time))
    }
    else {
        stop("The taxon name was not found in the tree data matrix. It should be a number between 1 and the number of nodes (internal and external)")
    }
  }



cleanEx()
nameEx("get.rate.descendant.pairs")
### * get.rate.descendant.pairs

flush(stderr()); flush(stdout())

### Name: get.rate.descendant.pairs
### Title: get.rate.descendant.pairs
### Aliases: get.rate.descendant.pairs
### Keywords: phylo

### ** Examples


set.seed(1234525)
myTree <- rcoal(50)

#Simulate rates with no autocorrelation (uncorrelated rates)

rateTreeUncor <- simulate.rate(myTree, FUN = simulate.uncor.lnorm)

uncorRates <- get.rate.descendant.pairs(rateTreeUncor)

#Simulate rates with high autocorrelation

rateTreeAutocor <- simulate.rate(myTree, FUN = simulate.autocor.kishino, params = list(initial.rate = 0.01, v = 0.01))

autocorRates <- get.rate.descendant.pairs(rateTreeAutocor)

# Plot the rates through time for all lineages to inspect the degree of autocorrelation

plot(rateTreeAutocor, col.lineages = rainbow(50))
plot(rateTreeUncor, col.lineages = rainbow(50))

# Estimate correlation of branch and daughter branch-wise rates. The correlation coefficient should be higher for the autocorrelated rates:

cor(x = autocorRates$parent.rate, y = autocorRates$daughter.rate)

cor(x = uncorRates$parent.rate, y = uncorRates$daughter.rate)


## The function is currently defined as
function (rate.sim.object) 
{
    dat <- rate.sim.object$tree.data.matrix
    parent.rate <- vector()
    daughter.rate <- vector()
    diff.br.len <- vector()
    for (i in 1:nrow(dat)) {
        daughter.temp <- dat[i, 3]
        parent.temp <- dat[i, 2]
        br.temp <- dat[i, 4]
        if (parent.temp %in% dat[, 3]) {
            parent.rate <- c(parent.rate, dat[dat[, 3] == parent.temp, 
                5])
            daughter.rate <- c(daughter.rate, dat[i, 5])
            diff.br.len <- c(diff.br.len, abs(br.temp - dat[dat[, 
                3] == parent.temp, 4]))
        }
    }
    return(data.frame(parent.rate, daughter.rate, diff.br.len))
  }



cleanEx()
nameEx("get.tree.data.matrix")
### * get.tree.data.matrix

flush(stderr()); flush(stdout())

### Name: get.tree.data.matrix
### Title: get.tree.data.matrix
### Aliases: get.tree.data.matrix
### Keywords: phylo

### ** Examples

set.seed(12345)
myTree <- rcoal(10)
myDataMatrix <- get.tree.data.matrix(myTree)
print(myDataMatrix)


## The function is currently defined as
function (phylo) 
{
    require(phangorn)
    require(geiger)
    data.matrix <- matrix(data = NA, ncol = 7, nrow = length(phylo$edge.length))
    colnames(data.matrix) <- c("branch.index", "parent.node", 
        "daughter.node", "branch.midage", "branch.rate", "length.subst", 
        "length.time")
    data.matrix[, 1] <- 1:length(phylo$edge.length)
    data.matrix[, 2] <- phylo$edge[, 1]
    data.matrix[, 3] <- phylo$edge[, 2]
    data.matrix[, 4] <- mid.edge.ages(phylo)
    data.matrix[, 7] <- phylo$edge.length
    class(data.matrix) <- "tree.data.matrix"
    return(data.matrix)
  }



cleanEx()
nameEx("mid.edge.ages")
### * mid.edge.ages

flush(stderr()); flush(stdout())

### Name: mid.edge.ages
### Title: mid.edge.ages obtans the ages of the mid point of the branches
###   of a phylogenetic tree. It is not necessary for the tree to be
###   ultrametric.
### Aliases: mid.edge.ages
### Keywords: phylo

### ** Examples

set.seed(12345)
myTree <- rcoal(10)

plot(myTree)
midAges <- mid.edge.ages(myTree)
edgelabels(round(midAges, 2)) 

## The function is currently defined as
function (phylo) 
{
    require(phangorn)
    rootage <- max(allnode.times(phylo))
    if (is.ultrametric(phylo) == TRUE) {
        midages <- vector()
        for (i in 1:length(phylo$edge.length)) {
            if (phylo$edge[i, 2] > length(phylo$tip.label)) {
                recent.node.age <- branching.times(phylo)[(phylo$edge[i, 
                  2] - length(phylo$tip.label))]
                halflength <- phylo$edge.length[i]/2
                midages[i] <- recent.node.age + halflength
            }
            else {
                midages[i] <- phylo$edge.length[i]/2
            }
        }
        return(midages)
    }
    else {
        nodetimes <- vector()
        extantedgelen <- max(phylo$edge.length[as.vector(which(phylo$edge[, 
            1] == as.numeric(names(which(branching.times(phylo) == 
            min(branching.times(phylo)))))))])
        addedval <- abs(min(branching.times(phylo))) + extantedgelen
        for (i in 1:length(branching.times(phylo))) {
            nodetimes[i] <- (rootage/(max(branching.times(phylo)) + 
                addedval)) * (branching.times(phylo) + addedval)[i]
        }
        brlen <- vector()
        for (i in 1:length(phylo$edge.length)) {
            brlen[i] <- (rootage/(max(branching.times(phylo)) + 
                addedval)) * phylo$edge.length[i]
        }
        midages <- vector()
        for (i in 1:length(brlen)) {
            if (phylo$edge[i, 2] > length(phylo$tip.label)) {
                daughter.node.age <- nodetimes[(phylo$edge[i, 
                  2] - length(phylo$tip.label))]
                halflength <- brlen[i]/2
                midages[i] <- daughter.node.age + halflength
            }
            else {
                parent.node.age <- nodetimes[(phylo$edge[i, 1] - 
                  length(phylo$tip.label))]
                midages[i] <- parent.node.age - (brlen[i]/2)
            }
        }
        return(round(midages, 5))
    }
  }



cleanEx()
nameEx("pathnode")
### * pathnode

flush(stderr()); flush(stdout())

### Name: pathnode
### Title: pathnode
### Aliases: pathnode
### Keywords: Node-density-effect

### ** Examples

set.seed(12345)
myTree <- rtree(10)

par(mfrow = c(1, 2))
plot(myTree)
nde <- pathnode(myTree)
nde



## The function is currently defined as
function (phylo, tipsonly = T) 
{
    require(phangorn)
    di.tr <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 
        2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, 
        ])
    if (tipsonly == TRUE) {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == 
            root.tr, 1:length(phylo$tip.label)]
        nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo, 
            x)))
    }
    else {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == 
            root.tr, ]
	nodesinpath <- sapply(1:(length(phylo$tip.label)+phylo$Nnode), function(x) length(Ancestors(phylo, x)))
    }
    plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", 
        ylab = "Number of parent nodes", pch = 20)
    return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
  }



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot.ratesim")
### * plot.ratesim

flush(stderr()); flush(stdout())

### Name: plot.ratesim
### Title: plot.ratesim
### Aliases: plot.ratesim
### Keywords: simulate.rate

### ** Examples


set.seed(123425)
myTree <- rcoal(10)
rateTree <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot.ratesim(rate.sim.object = rateTree, col.lineages = rainbow(10), type = "s")

# for a non-ultrametric tree
set.seed(123425)
myTree1 <- rtree(10)
rateTree1 <- simulate.rate(tree = myTree1, FUN = simulate.uncor.lnorm)
plot.ratesim(rate.sim.object = rateTree1, col.lineages = rainbow(10), type = "s")




## The function is currently defined as
function (rate.sim.object, col.lineages = colors(), type = "l") 
{
    rates.time.list <- list()
    for (i in 1:length(rate.sim.object[[1]]$tip.label)) {
        rates.time.list[[i]] <- get.lineage.time.rate(i, rate.sim.object)
    }
    ylims <- range(lapply(rates.time.list, function(y) range(y[, 
        2])))
    chrono <- rate.sim.object[[1]]
    chrono$edge.length <- rate.sim.object[[2]][, 7]
    node.ages <- allnode.times(chrono)
    xlims <- sort(range(node.ages), decreasing = T)
    par(mfrow = c(1, 2))
    plot(rates.time.list[[1]][, 1], rates.time.list[[1]][, 2], 
        ylim = ylims, xlim = xlims, ylab = "Rate", xlab = "Time", 
        type = type, lwd = 3, col = col.lineages[1])
    for (k in 2:length(rates.time.list)) {
        lines(rates.time.list[[k]][, 1], rates.time.list[[k]][, 
            2], ylim = ylims, xlim = xlims, col = col.lineages[k], 
            lwd = 3, type = type)
    }
    plot(chrono, edge.width = 1 + log(rate.sim.object[[2]][, 
        5]/min(rate.sim.object[[2]][, 5])), show.tip.label = F, 
        root.edge = T)
    tiplabels(pch = 16, col = col.lineages, cex = 1.5)
  }



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("simulate.autocor.kishino")
### * simulate.autocor.kishino

flush(stderr()); flush(stdout())

### Name: simulate.autocor.kishino
### Title: simulate.autocor.kishino
### Aliases: simulate.autocor.kishino
### Keywords: phylo rate of evolution

### ** Examples

set.seed(1234525)
myTree <- rcoal(20)

#Simulate high autocorrelation
kishinoRateTreeHigh <- simulate.autocor.kishino(myTree, params = list(initial.rate = 0.01, v = 0.001))
plot(kishinoRateTreeHigh, col.lineages = rainbow(20))

#Simulate low autocorrelation
kishinoRateTreeLow <- simulate.autocor.kishino(myTree, params = list(initial.rate = 0.01, v = 0.5))
plot(kishinoRateTreeLow, col.lineages = rainbow(20))

## The function is currently defined as
function (tree, params = list(initial.rate = 0.01, v = 0.3)) 
{
    require(phangorn)
    require(geiger)
    initial.rate <- params$initial.rate
    v = params$v
    data.matrix <- get.tree.data.matrix(tree)
    while (any(is.na(data.matrix[, 5])) | any(is.nan(data.matrix[, 
        5]))) {
        data.matrix[1, 5] <- initial.rate
        for (i in 2:nrow(data.matrix)) {
            parent.node <- data.matrix[i, 2]
            preceeding.parent <- data.matrix[, 2][data.matrix[, 
                3] == parent.node]
            preceeding.parent.brage <- data.matrix[, 4][data.matrix[, 
                2] == preceeding.parent][1]
            preceeding.parent.brrate <- data.matrix[, 5][data.matrix[, 
                2] == preceeding.parent][1]
            if (!(is.na(preceeding.parent.brrate)) & !(is.nan(preceeding.parent.brrate)) & 
                (parent.node %in% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), 
                  sd = v * data.matrix[i - 1, 7]^0.5))
            }
            else if (!(parent.node %in% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(initial.rate)), 
                  sd = sqrt(initial.rate)))
            }
        }
    }
    data.matrix[, 6] <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.autocor.thorne")
### * simulate.autocor.thorne

flush(stderr()); flush(stdout())

### Name: simulate.autocor.thorne
### Title: simulate.autocor.thorne
### Aliases: simulate.autocor.thorne
### Keywords: phylo rate of evolution

### ** Examples


set.seed(1234525)
myTree <- rcoal(20)

#Simulate high autocorrelation
thorneRateTreeHigh <- simulate.autocor.thorne(myTree, params = list(initial.rate = 0.01, v = 0.001))
plot(thorneRateTreeHigh, col.lineages = rainbow(20))

#Simulate low autocorrelation
thorneRateTreeLow <- simulate.autocor.thorne(myTree, params = list(initial.rate = 0.01, v = 0.5))
plot(thorneRateTreeLow, col.lineages = rainbow(20))



## The function is currently defined as
function (tree, params = list(initial.rate = 0.01, v = 0.3)) 
{
    require(phangorn)
    require(geiger)
    initial.rate <- params$initial.rate
    v = params$v
    data.matrix <- get.tree.data.matrix(tree)
    while (any(is.na(data.matrix[, 5])) | any(is.nan(data.matrix[, 
        5]))) {
        data.matrix[1, 5] <- initial.rate
        for (i in 2:nrow(data.matrix)) {
            parent.node <- data.matrix[i, 2]
            preceeding.parent <- data.matrix[, 2][data.matrix[, 
                3] == parent.node]
            preceeding.parent.brage <- data.matrix[, 4][data.matrix[, 
                2] == preceeding.parent][1]
            preceeding.parent.brrate <- data.matrix[, 5][data.matrix[, 
                2] == preceeding.parent][1]
            if (!(is.na(preceeding.parent.brrate)) & !(is.nan(preceeding.parent.brrate)) & 
                (parent.node %in% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), 
                  sd = v * abs(data.matrix[i, 4] - preceeding.parent.brage)^0.5))
            }
            else if (!(parent.node %in% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(initial.rate)), 
                  sd = sqrt(initial.rate)))
            }
        }
    }
    data.matrix[, 6] <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.clock")
### * simulate.clock

flush(stderr()); flush(stdout())

### Name: simulate.clock
### Title: simulate.clock
### Aliases: simulate.clock
### Keywords: clock phylo

### ** Examples

set.seed(1234525)

myTree <- rcoal(10)

# A tree with low stochastic variation
rateClock <- simulate.clock(tree = myTree, params = list(rate = 0.01, noise = 0.00001))

#Note the scale in the y axis. Rate variation is very low
plot(rateClock, col.lineages = rainbow(10))


## The function is currently defined as
function (tree, params = list(rate = 0.02, noise = 1e-04)) 
{
    rate <- params$rate
    noise <- params$noise
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rep(rate, times = length(tree$edge.length))
    branch.rates <- abs(branch.rates + rnorm(length(tree$edge.length), 
        mean = 0, sd = noise))
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.rate")
### * simulate.rate

flush(stderr()); flush(stdout())

### Name: simulate.rate
### Title: simulate.rate is the generic function to simulate rates along
###   phylogenetic trees with any of the models.
### Aliases: simulate.rate
### Keywords: phylo rate of evolution

### ** Examples

set.seed(1234525)

myTree <- rcoal(50)

# Simulate uncorrelated rates with default parameters:
rateTree.default <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot(rateTree.default, col.lineages = rainbow(50))

# Simulate uncorrelated rates with custom parameters:
rateTree.custom <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm, params = list(mean.log = -3.9, sd.log = 0.8))
plot(rateTree.custom, col.lineages = rainbow(50))


## The function is currently defined as
function (tree, FUN, ...) 
{
    ratesim.object <- FUN(tree, ...)
    return(ratesim.object)
  }



cleanEx()
nameEx("simulate.tdep.ho")
### * simulate.tdep.ho

flush(stderr()); flush(stdout())

### Name: simulate.tdep.ho
### Title: simulate.tdep.ho
### Aliases: simulate.tdep.ho
### Keywords: phylogenetics molecular rates time-dependence

### ** Examples


set.seed(12345)
myTree <- rcoal(50)
plot(myTree); axisPhylo()
# Perhaps a lamda value of 4 is more appropriate to simulate a significant rate change through time in this phylogeny. We also add additional noise, and leave the other default values unchanged.
plot(function(x) 0.015 + (0.035 * exp(-4 * x)), xlim = c(0, max(branching.times(myTree))))
rate.simulation <- simulate.tdep.ho(myTree, params = list(mu = 0.035, srate = 0.015, lambda = 4, noise = 0.003))
plot(rate.simulation[[2]][,4], rate.simulation[[2]][,5], pch = 19, xlab = "Mid age of branch", ylab = "Molecular rate")
plot(rate.simulation[[1]])

## The function is currently defined as
function (tree, params = list(mu = 0.035, srate = 0.015, lambda = 0.1, 
    noise = 0.001)) 
{
    require(phangorn)
    require(geiger)
    mu <- params$mu
    srate <- params$srate
    lambda <- params$lambda
    noise <- params$noise
    fun.rate <- function(x, m = mu, s = srate, lam = lambda) {
        if (any(x >= 0)) {
            return(s + (m * exp(-lam * x)))
        }
        else {
            stop("x is cannot be a negative number")
        }
    }
    data.matrix <- get.tree.data.matrix(tree)
    node.ages <- allnode.times(tree)
    b.times <- c(rep(0, length(tree$tip.label)), node.ages[(length(tree$tip.label) + 
        1):length(node.ages)])
    names(b.times) <- 1:length(b.times)
    ratetemp <- vector()
    for (i in 1:length(tree$edge.length)) {
        parentage <- b.times[as.character(data.matrix[i, 2])]
        daughterage <- b.times[as.character(data.matrix[i, 3])]
        ratetemp[i] <- integrate(fun.rate, lower = daughterage, 
            upper = parentage)$value/data.matrix[i, 7]
    }
    data.matrix[, 5] <- abs(ratetemp + rnorm(nrow(data.matrix), 
        mean = 0, sd = noise))
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.uncor.exp")
### * simulate.uncor.exp

flush(stderr()); flush(stdout())

### Name: simulate.uncor.exp
### Title: simulate.uncor.exp
### Aliases: simulate.uncor.exp
### Keywords: rate of evolution phylo

### ** Examples

set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.uncor.exp(tree = myTree, params = list(mean.exp = 0.01))
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])

## The function is currently defined as
function (tree, params = list(mean.exp = 0.001)) 
{
    mean.exp <- params$mean.exp
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rexp(n = length(tree$edge.length), rate = 1/mean.exp)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.uncor.gamma")
### * simulate.uncor.gamma

flush(stderr()); flush(stdout())

### Name: simulate.uncor.gamma
### Title: simulate.uncor.gamma
### Aliases: simulate.uncor.gamma
### Keywords: rate of evolution phylo

### ** Examples


set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.uncor.gamma(tree = myTree, params = list(shape = 98, rate = 4361))
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])

## The function is currently defined as
function (tree, params = list(shape = 98, rate = 4361)) 
{
    shape.gamma <- params$shape
    rate.gamma <- params$rate
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rgamma(n = length(tree$edge.length), shape = shape.gamma, 
        rate = rate.gamma)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.uncor.lnorm")
### * simulate.uncor.lnorm

flush(stderr()); flush(stdout())

### Name: simulate.uncor.lnorm
### Title: simulate.uncor.lnorm
### Aliases: simulate.uncor.lnorm
### Keywords: rate of evolution phylo

### ** Examples

set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.uncor.lnorm(tree = myTree, params = list(mean.log = -3.9, sd.log = 0.5))
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])


## The function is currently defined as
function (tree, params = list(mean.log = -3.9, sd.log = 0.1)) 
{
    mean.log <- params$mean.log
    sd.log <- params$sd.log
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rlnorm(n = length(tree$edge.length), meanlog = mean.log, 
        sdlog = sd.log)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("simulate.white.noise")
### * simulate.white.noise

flush(stderr()); flush(stdout())

### Name: simulate.white.noise
### Title: simulate.white.noise
### Aliases: simulate.white.noise
### Keywords: white noise phylo

### ** Examples



set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.white.noise(tree = myTree, params = list(mean.log = -3.9, sd.log = 0.1))
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])

## The function is currently defined as
function (tree, params = list(mean.log = -3.9, sd.log = 0.1)) 
{
    mean.log <- params$mean.log
    sd.log <- params$sd.log
    data.matrix <- get.tree.data.matrix(tree)
    branch.noise <- sd.log / (data.matrix[, 7]/mean(data.matrix[, 
        7]))
    branch.rates <- sapply(branch.noise, function(x) rlnorm(1, 
        mean.log, branch.noise))
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, tree.data.matrix = data.matrix)
    class(res) <- "ratesim"
    return(res)
  }



cleanEx()
nameEx("trann2trdat")
### * trann2trdat

flush(stderr()); flush(stdout())

### Name: trann2trdat
### Title: trann2trdat. extract annontations from a nexus tree.
### Aliases: trann2trdat
### Keywords: Phylo

### ** Examples

## Not run: 
##D myAnnotatedTree <- read.annotated.nexus("annotated.tree")
##D annnotationData <- trann2trdat(myAnnotatedTree)
##D head(annotationData)
## End(Not run)
## The function is currently defined as
function (tree) 
{
    require(epibase)
    tree$edge.length <- unlist(sapply(tree$annotations, function(x) {
        x$length
    }))[1:length(tree$edge.length)]
    rates <- unlist(sapply(tree$annotations, function(x) {
        x$rate_median
    }))
    if (is.ultrametric(tree) == TRUE) {
        midages <- mid.edge.ages(tree)
    }
    else {
        midages <- mid.edge.ages(tree, max(unlist(sapply(tree$annotations, 
            function(x) {
                x$height_median
            }))))
    }
    timelen <- tree$edge.length
    subslen <- tree$edge.length * rates
    return(data.frame(branch = rownames(as.data.frame(tree$edge)), 
        parent = tree$edge[, 1], daughter = tree$edge[, 2], midage = midages, 
        rate = rates, blensubs = subslen, blentime = timelen))
  }



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
