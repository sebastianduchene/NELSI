\name{get.lineage.time.rate}
\alias{get.lineage.time.rate}
\title{
get.lineage.time.rate
}
\description{
get.lineage.time.rate obtains the branch-wise rates along a lineage. It uses a taxon to trace the rates back to the root of the tree.
It is used internally by plot.ratesim().
}
\usage{
get.lineage.time.rate(taxon, sim.rate.object)
}
\arguments{
  \item{taxon}{
An integer corresponding to the index of a taxon. The index is the order of the taxa in the tree, which can be obtained with tr$tip.label. See example.
}
  \item{sim.rate.object}{
An object of class ratesim. Usually produced with simulate.rate
}
}
\details{
None
}
\value{
A data frame with two colums. The rows correspond to the each of the branches from the taxon to the root of the tree.
The first column of the data frame is the mid age of the branch, and the second column is the rate along the branch.
}
\references{
None.
}
\author{
Sebastian Duchene
}
\note{
None
}

\seealso{
simulate.rate()
}
\examples{
set.seed(123425)
myTree <- rcoal(10)
rateTree <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot(rateTree, col.lineages = rainbow(10))

# Get the rate for the lineages ascending from the first taxon
# Find the name of the first taxon
myTree$tip.label

get.lineages.rate(taxon = 1, sim.rate.object = rateTree)



## The function is currently defined as
function (taxon, sim.rate.object) 
{
    tree.data.matrix <- sim.rate.object[[2]]
    chrono <- sim.rate.object[[1]]
    chrono$edge.length <- tree.data.matrix[, 7]
    taxon.init <- taxon
    if (taxon \%in\% tree.data.matrix[, 3]) {
        data.matrix <- tree.data.matrix
        branch.times <- vector()
        rate.time <- vector()
        repeat {
            parent <- data.matrix[, 2][data.matrix[, 3] == taxon]
            time.br <- data.matrix[, 4][data.matrix[, 3] == taxon]
            rate.br <- data.matrix[, 5][data.matrix[, 3] == taxon]
            rate.time <- c(rate.time, rate.br)
            branch.times <- c(branch.times, time.br)
            taxon <- parent
            if (!(parent \%in\% data.matrix[, 3])) {
                break
            }
        }
        first.rate <- rate.time[length(rate.time)]
        last.rate <- rate.time[1]
        rate.time <- c(last.rate, rate.time, first.rate)
        node.times <- allnode.times(chrono)
        root.age <- max(node.times)
        branch.times <- c(node.times[taxon.init], branch.times, 
            root.age)
        return(data.frame(branch.times, rate.time))
    }
    else {
        stop("The taxon name was not found in the tree data matrix. It should be a number between 1 and the number of nodes (internal and external)")
    }
  }
}
\keyword{ substitution rate }
\keyword{ phylo }
