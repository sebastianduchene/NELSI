\name{simulate.autocor.thorne}
\alias{simulate.autocor.thorne}
\title{
simulate.autocor.thorne
}
\description{
Simulate rates of evolution along a phylogenetic tree with the model reported in Thorne et al. (1998).
}
\usage{
simulate.autocor.thorne(tree, params = list(initial.rate = 0.01, v = 0.3))
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
Thorne, J.L., Kishino,H., and Painter, I.S. "Estimating the rate of evolution of the rate of molecular evolution." Molecular Biology and Evolution 15.12 (1998): 1647-1657.
}
\author{
Sebastian Duchene. See the reference for the method.
}
\note{
None
}

\seealso{
simulate.autocor.kishino for a similar model
}
\examples{

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
                (parent.node \%in\% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), 
                  sd = v * abs(data.matrix[i, 4] - preceeding.parent.brage)^0.5))
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
