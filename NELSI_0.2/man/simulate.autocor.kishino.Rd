\name{simulate.autocor.kishino}
\alias{simulate.autocor.kishino}
\title{
simulate.autocor.kishino
}
\description{
Simulate rates of evolution along a phylogenetic tree with the model reported in Kishino et al. (2001).
}
\usage{
simulate.autocor.kishino(tree, params = list(initial.rate = 0.01, v = 0.3))
}
\arguments{
  \item{tree}{
A phylogeneti tree of class 'phylo'. The tree should be a chronogram, with branch lengths in time.
}
  \item{params}{
parameters for the autocorrelation function. This should be a list with two items:
 \item{inital.rate}{The rate at the root of the tree}
 \item{v}{This is the nu parameter described in Kishino et al. (2001). A high value implies low autocorrelation, with high differences in the rates of parent and daughter branches. A low value of v implies high autocorrelation, with very similar rates between parent and daughter branches.
}
}
}

\details{
See the original reference for further details.
}
\value{
An object of class 'ratesim'. This is a list with two items:
 \item{phylogram }{The phylogenetic tree with branch lengths in units of substitutions (phylogram)}
 \item{tree.data.matrix }{This is a matrix with the number of substitutions, rates, and times along every branch in the tree. See get.tree.data.matrix for more details}
}
\references{
Kishino, H., Thorne, J.L., and Bruno, W. J. "Performance of a divergence time estimation method under a probabilistic model of rate evolution." Molecular Biology and Evolution 18.3 (2001): 352-361.
}
\author{
Sebastian Duchene. But see the reference for the original method.
}
\note{
None
}

\seealso{
simulate.rate() can call all the rate simulation functions internally. 
}
\examples{
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
                (parent.node \%in\% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), 
                  sd = v * data.matrix[i - 1, 7]^0.5))
            }
            else if (!(parent.node \%in\% data.matrix[, 3])) {
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
}
\keyword{ phylo }
\keyword{ rate of evolution }
