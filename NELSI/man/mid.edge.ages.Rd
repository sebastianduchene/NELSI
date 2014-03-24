\name{mid.edge.ages}
\alias{mid.edge.ages}
\title{
mid.edge.ages obtans the ages of the mid point of the branches of a phylogenetic tree. It is not necessary for the tree to be ultrametric. 
}
\description{
mid.edge.ages produces a vector of length = number of edges of the tree with the mid.age of the branchies. Assuming that the tree is a chronogram.
}
\usage{
mid.edge.ages(phylo)
}
\arguments{
  \item{phylo}{
A phylogenetic tree of class 'phylo'. The branch lengths should be units of time.
}
}
\details{
None
}
\value{
A vector of length = number of edges in the tree. Each value is the median branch age.
}
\references{
None
}
\author{
David Duchene
}
\note{
%%  ~~further notes~~
None
}


\seealso{
allnode.times produces the ages of internal nodes and tips.
}
\examples{
set.seed(12345)
myTree <- rcoal(10)

plot(myTree)
midAges <- mid.edge.ages(myTree)
edgelabels(round(midAges, 2)) 

## The function is currently defined as
function (phylo) 
{
    require(phangorn)
    rootage <- max(allnode.times(phylo))
    if (is.ultrametric(phylo) == TRUE) {
        midages <- vector()
        for (i in 1:length(phylo$edge.length)) {
            if (phylo$edge[i, 2] > length(phylo$tip.label)) {
                recent.node.age <- branching.times(phylo)[(phylo$edge[i, 
                  2] - length(phylo$tip.label))]
                halflength <- phylo$edge.length[i]/2
                midages[i] <- recent.node.age + halflength
            }
            else {
                midages[i] <- phylo$edge.length[i]/2
            }
        }
        return(midages)
    }
    else {
        nodetimes <- vector()
        extantedgelen <- max(phylo$edge.length[as.vector(which(phylo$edge[, 
            1] == as.numeric(names(which(branching.times(phylo) == 
            min(branching.times(phylo)))))))])
        addedval <- abs(min(branching.times(phylo))) + extantedgelen
        for (i in 1:length(branching.times(phylo))) {
            nodetimes[i] <- (rootage/(max(branching.times(phylo)) + 
                addedval)) * (branching.times(phylo) + addedval)[i]
        }
        brlen <- vector()
        for (i in 1:length(phylo$edge.length)) {
            brlen[i] <- (rootage/(max(branching.times(phylo)) + 
                addedval)) * phylo$edge.length[i]
        }
        midages <- vector()
        for (i in 1:length(brlen)) {
            if (phylo$edge[i, 2] > length(phylo$tip.label)) {
                daughter.node.age <- nodetimes[(phylo$edge[i, 
                  2] - length(phylo$tip.label))]
                halflength <- brlen[i]/2
                midages[i] <- daughter.node.age + halflength
            }
            else {
                parent.node.age <- nodetimes[(phylo$edge[i, 1] - 
                  length(phylo$tip.label))]
                midages[i] <- parent.node.age - (brlen[i]/2)
            }
        }
        return(round(midages, 5))
    }
  }
}

\keyword{ phylo }
