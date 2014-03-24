\name{simulate.uncor.lnorm}
\alias{simulate.uncor.lnorm}
\title{
simulate.uncor.lnorm
}
\description{
Simulate the rate of evolution along a phylogenetic tree using the lognormal model described in Drummond et al. (2006)
}
\usage{
simulate.uncor.lnorm(tree, params = list(mean.log = -3.9, sd.log = 0.1))
}
\arguments{
  \item{tree}{A phylogenetic tree of class 'phylo'. The branch lengths should be in units of time (chronogram)}
\item(params}{
parameters for the uncorrelated rates function. This should be a list with two items:
 \item{mean.log}{This is the mean rate of the lognormal distribution. It should be in log scale}
 \item{sd.log}{The standard deviation of the lognormal distribution}
}
}
\details{
None
}
\value{
An object of class 'ratesim'. This is a list with two items:
 \item{phylogram }{The phylogenetic tree with branch lengths in units of substitutions (phylogram)}
 \item{tree.data.matrix }{This is a matrix with the number of substitutions, rates, and times along every branch in the tree. See get.tree.data.matrix for more details}
}
}
\references{
Drummond, A.J., et al. "Relaxed phylogenetics and dating with confidence." PLoS biology 4.5 (2006): e88.
}
\author{
Sebastian Duchene. See the reference for the description of the model.
}
\note{
None
}


\seealso{
simulte.uncor.exp
}
\examples{
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
}
\keyword{ rate of evolution }
\keyword{ phylo }
