\name{simulate.tdep.generic}
\alias{simulate.tdep.generic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate.tdep.generic
}
\description{
This function simulates time-dependent molecular rates for branches in a chronogram using any specified function of rate through time.
}
\usage{
simulate.tdep.generic(tree, noise = 0.001, fun = function(time) {
    0.015 * (0.035 * exp(-0.1 * time))
})
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree of class 'phylo', with branch lengths in terms of time.
}
  \item{noise}{
This is the standard deviation of a normal distribution from which noise for the rate is added.
}
  \item{fun}{
A function of the rate through time. The default is a vertically transformed exponential function. All the parameters for the function must be provided, except for time.
}
}
\details{
This function generalises simulate.tdep.ho so that the user can provide any function for the (time-dependent) variation of the rate through time.
}
\value{
An object of class "ratesim". This is similar to a list; the first element is of class "phylo" with branch lengths in terms of substitutions. The second is a "tree.data.matrix" which can be used as a "data.frame". The tree.data.matrix contains the edge object of the phylogeny, the mid age of a branch, and branch lengths in terms of substitutions, time, and molecular rate.
}
\references{
Pending.
}
\author{
David Duchene and Sebastian Duchene
}
\note{
None
}
\seealso{
simulate.tdep.ho
simulate.autocor.kishino
simulate.autocor.thorne
simulate.uncor.exp
simulate.uncor.gamma
simulate.uncor.lnorm
simulate.white.noise
simulate.clock
}
\examples{

set.seed(12345)
myTree <- coal(50)
plot(myTree); axisPhylo()
# Using a simple linear function of rate variation:
rate.simulation <- simulate.tdep.generic(myTree, noise = 0.0001, fun = function(x) -0.015*x + 0.035)
plot(rate.simulation[[2]][,4], rate.simulation[[2]][,5], pch = 19, xlab = "Mid age of branch", ylab = "Molecular rate")
plot(rate.simulation[[1]])

## The function is currently defined as
function (tree, noise = 0.001, fun = function(time) {
    0.035 * (0.015 * exp(-0.1 * time))
  }) 
{
    require(phangorn)
    require(geiger)
    fun.rate <- function(time) {
        if (any(time >= 0)) {
            return(fun(time))
        }
        else {
            stop("x is cannot be a negative number")
        }
    }
    data.matrix <- get.tree.data.matrix(tree)
    node.ages <- allnode.times(tree)
    b.times <- allnode.times(tree)
    b.times <- c(rep(0, length(tree$tip.label)), node.ages[(length(tree$tip.labels) + 
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
}

\keyword{ phylogenetics }
\keyword{ molecular rates }
\keyword{ time-dependence }

