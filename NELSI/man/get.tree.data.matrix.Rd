\name{get.tree.data.matrix}
\alias{get.tree.data.matrix}
\title{
get.tree.data.matrix
}
\description{
Obtain the tree data matrix from a phylogenetic tree. The tree should be a chronogram, with branch lengths representing time.
}
\usage{
get.tree.data.matrix(phylo)
}
\arguments{
  \item{phylo}{
A phylogenetic tree of class 'phylo'. The tree should be a chronogram, with branch lengths representing time. 
}
}
\details{
None
}
\value{
A matrix with 7 colums and number of rows = number of edges in the tree. Each row corresponds to a branch (or edge) of the tree. The columns of the matrix correspond to the following:

branch.index: The number of the branch (or edge) in the tree.
parent.node: The node of the corresponding branch that is closer to the root of the tree.
daughter.node: The node of the corresponding branch that is closer to the tips.
branch.midage: The median age of the branch.
branch.rate: The rate along the branch. It is set to NA if get.tree.data.matrix is called directly on a phylogenetic trees.  The functions simulate.rate and trann2trdat fill this column with the branch-wise rates.
lengths.subst: Number of substitutions along the branch. It is set to NA if get.tree.data.matrix is called directly on a phylogenetic trees.  The functions simulate.rate and trann2trdat fill this column with the branch-wise rates.
length.time: The branch lengths in units of time.
}
\references{
None.
}
\author{
David Duchene and Sebastian Duchene
}
\note{
This function is used internally by the simulate.rate()]
}

\seealso{
simulate.rate
}
\examples{
set.seed(12345)
myTree <- rcoal(10)
myDataMatrix <- get.tree.data.matrix(myTree)
print(myDataMatrix)


## The function is currently defined as
function (phylo) 
{
    require(phangorn)
    require(geiger)
    data.matrix <- matrix(data = NA, ncol = 7, nrow = length(phylo$edge.length))
    colnames(data.matrix) <- c("branch.index", "parent.node", 
        "daughter.node", "branch.midage", "branch.rate", "length.subst", 
        "length.time")
    data.matrix[, 1] <- 1:length(phylo$edge.length)
    data.matrix[, 2] <- phylo$edge[, 1]
    data.matrix[, 3] <- phylo$edge[, 2]
    data.matrix[, 4] <- mid.edge.ages(phylo)
    data.matrix[, 7] <- phylo$edge.length
    class(data.matrix) <- "tree.data.matrix"
    return(data.matrix)
  }
}
\keyword{ phylo }
