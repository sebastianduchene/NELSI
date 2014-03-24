\name{simulate.clock}
\alias{simulate.clock}
\title{
simulate.clock
}
\description{
simulate.clock simulates a constant rate of evolution along a phylogenetic tree
}
\usage{
simulate.clock(tree, params = list(rate = 0.02, noise = 1e-04))
}
\arguments{
  \item{tree}{
A phylogeneti tree of class 'phylo'. The tree should be a chronogram, with branch lengths in time.
}
  \item{params}{
parameters for the clock rate simulation function. This should be a list with two items:
 \item{rate}{The rate for the tree}
 \item{noise}{The stochastic noise. Note that if this parameter is set very high compared to the rate, the model will be similar to an uncorrelated rates model}
}
\details{
None
}
\value{
An object of class 'ratesim'. This is a list with two items:
 \item{phylogram }{The phylogenetic tree with branch lengths in units of substitutions (phylogram)}
 \item{tree.data.matrix }{This is a matrix with the number of substitutions, rates, and times along every branch in the tree. See get.tree.data.matrix for more details}
}
\references{
This is a constant molecular rate model. The original model can be found in:
Zuckerkandl, Emile, and Linus Pauling. "Molecular disease, evolution and genetic heterogeneity." (1962): 189-225.

Zuckerkandl, Emile, and Linus Pauling. "Evolutionary divergence and convergence in proteins." Evolving genes and proteins 97 (1965): 97-166.
}
\author{
Sebastian Duchene. See references for the original description of the molecular clock.
}
\note{
None
}

\seealso{
None.
}
\examples{
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clock }
\keyword{ phylo }

