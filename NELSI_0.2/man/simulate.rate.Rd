\name{simulate.rate}
\alias{simulate.rate}
\title{
simulate.rate is the generic function to simulate rates along phylogenetic trees with any of the models.
}
\description{
simulate.rate simulates the rate with any of the models implemented in NELSI
}
\usage{
simulate.rate(tree, FUN, ...)
}
\arguments{
  \item{tree}{
A phylogenetic tree of class 'phylo'. The branch lengths should be in units of time (a chronogram)}
  \item{FUN}{
This is any of the rate simulation functions (please see the help for each for details on the models):
 - simulate.autocor.kishino
 - simulate.autocor.thorne
 - simulate.clock
 - simulate.tdep.generic
 - simulate.tdep.ho
 - simulate.uncor.exp
 - simulate.uncor.lnorm
 - simulate.uncor.gamma
 - simulate.white.noise
}
  \item{\dots}{
This should be the parameters for the rate models. To specify this use: params = list(parameter1 , parameter2). See the help files for each rate simulation function for details on the parameters.
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
\references{
See the original reference for NELSI: Pending.
}
\author{
David Duchene and Sebastian Duchene
}
\note{
None.
}

\seealso{
 - simulate.autocor.kishino
 - simulate.autocor.thorne
 - simulate.clock
 - simulate.tdep.generic
 - simulate.tdep.ho
 - simulate.uncor.exp
 - simulate.uncor.lnorm
 - simulate.uncor.gamma
 - simulate.white.noise
}
\examples{
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
}
\keyword{ phylo }
\keyword{ rate of evolution }