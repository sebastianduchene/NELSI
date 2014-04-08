\name{simulate.white.noise}
\alias{simulate.white.noise}
\title{
simulate.white.noise
}
\description{
simulate.white.noise simulates the rate along a phylogeneti tree according to a white noise process. In particular, it samples from a lognormal distribution with a given mean and standard deviation. For each branch a rate is drawn from a distribution with the given mean. The standard deviation is divided by the branch length divided by the mean of the branch lengths. As a result, the rate of a long branch will be drawn from a distribution with a lower standard deviation than that for a short branch (i.e. the rate is more uncertain for shorter branches). 
}
\usage{
simulate.white.noise(tree, params = NULL)
}
\arguments{
  \item{tree}{
A phylogenetic tree of class 'phylo'. The branch lengths should be in units of time (a chronogram).
}
  \item{params}{A list with the parameters for the simulation function, corresponding to the following two items:}
  \item{mean.log}{The mean of the lognormal distribution. It should be in log scale}
  \item{sd.log}{The standard deviation of the lognormal distribution. The actual value used for each branch will vary inversely with the branch length. See the description of the function}
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
A similar model is described in:
Lepage, Thomas, et al. "A general comparison of relaxed molecular clock models." Molecular biology and evolution 24.12 (2007): 2669-2680.
}
\author{
Sebastian Duchene and David Duchene
}
\note{
none
}


\seealso{
simulate.uncor.lnorm
}
\examples{


set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.white.noise(tree = myTree, params = NULL)
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])

## The function is currently defined as
function (tree, params = NULL) 
{
    data.matrix <- get.tree.data.matrix(tree)
    means.rates <- data.matrix[, 7] / sum(data.matrix[, 7])
    branch.rates <- sapply(means.rates, function(x) rlnorm(1, log(x), x))
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, tree.data.matrix = data.matrix)
    class(res) <- "ratesim"
    return(res)
  }
}
\keyword{ white noise }
\keyword{ phylo }
