\name{allnode.times}
\alias{allnode.times}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
allnode.times
}
\description{
allnode.times is used to obtain the ages of all nodes in the tree. It is similar to branching.times() in ape, but it also returns the ages of the tips. This is particularly useful for heterochronous trees.
}
\usage{
allnode.times(phylo, tipsonly = FALSE)
}
\arguments{
  \item{phylo}{
A phylogenetic tree of class 'phylo'
}
  \item{tipsonly}{
A logical (T / F). With T, the function returns the ages of the tips only. With F, it returns the ages of the tips and internal nodes.
}
}
\details{
This function is similar to branching.times, but it also returns the ages of the tips. In ultrametric trees, the ages of all tips are 0, but in heterocrhonous trees they can have different ages.
}
\value{
A vector with the ages of all nodes in the tree. The youngest tip always has age 0.
The items of the vector are numbered according to the 'phylo' object. To inspect the numbering of the tips and internal nodes see the example.
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
branching.times() from package 'ape'
}
\examples{
set.seed(12345)
myTree <- rtree(10)
plot(myTree)
# See the numbering of internal nodes and tips. Note that the tip labels and the actual numbering of the tips are different.
nodelabels()
tiplabels()
allnode.times(myTree)

# Plot the tree and add the ages of the tips and internal nodes.
plot(myTree, show.tip.label = F)
allTimes <- allnode.times(myTree)
allTimes
tiplabels(round(allTimes[1:10]))
nodelabels(round(allTimes[11:19]))


## The function is currently defined as
function (phylo, tipsonly = FALSE) 
{
    require(phangorn)
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] \%in\% phylo$edge[, 
        2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == 
        root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == 
        root.phylo, ]
    if (tipsonly == TRUE) {
        node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == 
            root.phylo, 1:length(phylo$tip.label)]
    }
    return(node.times)
  }
}
\keyword{ chronogram }
