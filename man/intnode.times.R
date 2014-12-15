\name{intnode.times}
\alias{itno.times}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
allnode.times
}
\description{
intnode.times is used to obtain the ages of all nodes in the tree. It is similar to branching.times() in ape, but it also returns the ages of the tips. This is particularly useful for heterochronous trees.
}
\usage{
intnode.times(phylo, tipsonly = FALSE)
}
\arguments{
  \item{phylo}{
A phylogenetic tree of class 'phylo'
}
}
\details{
This function is similar to branching.times, but it supports non ultrametric trees. The youngest tip always has an age of 0.
}
\value{
A vector with the ages of all internal nodes in the tree.
}
\references{
Pending.
}
\author{
David Duchene and Sebastian Duchene
}
\note{
None
}

\seealso{
branching.times() from package 'ape', and allnode.times from 'NELSI'.
}
\examples{
set.seed(12345)
myTree <- rtree(10)
plot(myTree)
# See the numbering of internal nodes and tips. Note that the tip labels and the actual numbering of the tips are different.
nodelabels()
tiplabels()
intnode.times(myTree)

# Plot the tree and add the ages of the tips and internal nodes.
plot(myTree, show.tip.label = FALSE)
nodetimes <- allnode.times(myTree)
nodetimes
tiplabels(round(nodetimes[1:10]))
nodelabels(round(nodetimes[11:19]))


## The function is currently defined as

intnode.times <-
function(phylo, tipsonly = FALSE){
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, (length(phylo$tip.label) + 1):nrow(di.phylo)]
    return(node.times)
}


}
\keyword{ chronogram }
