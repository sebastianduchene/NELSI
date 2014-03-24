\name{simulate.tdep.ho}
\alias{simulate.tdep.ho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate.tdep.ho
}
\description{
This function simulates time-dependent molecular rates for branches in a chronogram using an vertically transformed exponential curve (Ho et al. 2005).
}
\usage{
simulate.tdep.ho(tree, params = list(mu = 0.015, srate = 0.035, lambda = 0.1, noise = 0.001))
}
\arguments{
  \item{tree}{
A phylogenetic tree of class 'phylo', with branch lengths in terms of time.
}
  \item{params}{
These are the three parameters describing the variation in rate through time. mu is the instantaneous mutation rate, observed near the present. srate is the long term observed substitutions rate. lambda is the rate of the rate decrease towards the past. noise is the standard deviation of a normal distribution from which noise for the rate is added.
}
}
\details{
The number of substitutions for each branch is calculated by integrating the function given by the user within the limits of the age of a branch. Values of substitutions are added noise and used to determine the rate at each branch.

Note that this is function only uses one of the possible models for time dependence (see simulate.tdep.generic for a more general form).
}
\value{
An object of class "ratesim". This is similar to a list; the first element is of class "phylo" with branch lengths in terms of substitutions. The second is a "tree.data.matrix" which can be used as a "data.frame". The tree.data.matrix contains the edge object of the phylogeny, the mid age of a branch, and branch lengths in terms of substitutions, time, and molecular rate.
}
\references{
Ho, S. Y., Phillips, M. J., Cooper, A., & Drummond, A. J. (2005). Time dependency of molecular rate estimates and systematic overestimation of recent divergence times. Molecular biology and evolution, 22(7), 1561-1568.
}
\author{
David Duchene and Sebastian Duchene
}
\note{
None
}
\seealso{
simulate.tdep.generic
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
}

\keyword{ phylogenetics }
\keyword{ molecular rates }
\keyword{ time-dependence }
