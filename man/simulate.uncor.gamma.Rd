\name{simulate.uncor.gamma}
\alias{simulate.uncor.gamma}
\title{
simulate.uncor.gamma
}
\description{
Simulate the rate of evolution along a phylogenetic tree using the gamma distributed rates model described in Drummond et al. (2006)
}
\usage{
simulate.uncor.gamma(tree, params = list(shape = 3.98, rate = 516.53))
}
\arguments{

  \item{tree}{A phylogenetic tree of class 'phylo'. The branch lengths should be in units of time (chronogram)}
\item{params}{
parameters for the uncorrelated rates function. This should be a list with two items:}
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
Drummond, A.J., et al. "Relaxed phylogenetics and dating with confidence." PLoS biology 4.5 (2006): e88.
}
\author{
Sebastian Duchene. See the reference for the original description of the model.
}
\note{
Notes
}

\seealso{
simulate.uncor.lnorm
simulate.uncor.exp
}
\examples{

set.seed(1234525)

myTree <- rcoal(50)

rateTree <- simulate.uncor.gamma(tree = myTree, params = list(shape = 3.98, rate = 516.53))
plot(rateTree, col.lineages = rainbow(50))

#See the histogram of the branch-wise rates
hist(rateTree$tree.data.matrix[, 5])

## The function is currently defined as
function (tree, params = list(shape = 3.98, rate = 516.53)) 
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
}
\keyword{ rate of evolution }
\keyword{ phylo }
