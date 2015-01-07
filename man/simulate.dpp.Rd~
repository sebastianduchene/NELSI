\name{simulate.dpp}
\alias{simulate.dpp}
\title{
simulate.dpp
}
\description{
Simulate the rate of evolution along a phylogenetic tree using the dirichlet process described by Heath et al. (2012).
}
\usage{
simulate.dpp(tree, params = list(alpha = 1, shape = 10, rate = 0.0001))
}
\arguments{

  \item{tree}{A phylogenetic tree of class 'phylo'. The branch lengths should be in units of time (chronogram)}
\item{params}{
parameters for the rates function. This should be a list with two items:}
 \item{alpha}{The concentration of the dirichlet distribution. Lower values will cause more rate categories and rate heterogeneity.}
 \item{shape}{The shape of the gamma distribution. It follows the usual parametrisation of the gamma distribution}
 \item{rate}{The rate of the gamma distribution. Note that this is different from the substitution rate. Instead, this the the rate parameter typically used for the gamma distribution}
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
Heath, T. A., Holder, M. T., & Huelsenbeck, J. P. (2012). A dirichlet process prior for estimating lineage-specific substitution rates. Molecular biology and evolution, 29(3), 939-955.
}
\author{
David Duchene. See the reference for the original description of the model.
}
\note{
Notes
}

\seealso{
simulate.uncor.gamma
simulate.uncor.exp
}
\examples{

set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.dpp(tree = myTree, params = list(alpha = 1, shape = 10, rate = 0.0001))
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])

## The function is currently defined as

simulate.dpp <-
function(tree, params = list(alpha = 1, shape = 10, rate = 0.0001)){
    shape.gamma <- params$shape
    rate.gamma <- params$rate
    alpha <- params$alpha
    data.matrix <- get.tree.data.matrix(tree)
    nbranches <- nrow(tree$edge)
    cats <- sample(1:nbranches, 1, prob = rdirichlet(1, rep(alpha, nbranches)))
    branch.cats <- sample(1:cats, nbranches, replace = T, prob = rdirichlet(1, rep(alpha, cats)))
    branch.rates <- rgamma(n = cats, shape = shape.gamma, rate = rate.gamma)
    names(branch.rates) <- 1:cats
    branch.rates <- branch.rates[branch.cats]
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}

}
\keyword{ rate of evolution }
\keyword{ phylo }