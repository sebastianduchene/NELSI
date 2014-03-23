\name{pathnode}
\alias{pathnode}
\title{
pathnode
}
\description{
pathnode obtains the root-to tip distance for all tips in the tree, and the number of nodes along the path to each tip. These data are useful to evaluate the Node Density effect in phylogenetic trees.
}
\usage{
pathnode(phylo, tipsonly = T)
}
\arguments{
  \item{phylo}{
A phylogenetic tree of class 'phylo'
}
  \item{tipsonly}{
A logical value (T / F) indicating whether the function should only return the root-to-tip distance and number of nodes along a lineage for the tips only. Select T for the tips only, and F for the tips and internal nodes.
}
}
\details{
The function plots the root-to-tip distance and the number of nodes along each lineage. It also returns a list with the two variables. The data can be used to assess the node-density effect in phylogenetic trees.
}
\value{
A list with two items:
  \item{roottotippath }{A vector with the distance from the root to the tips. If tipsonly==F, it will also include the distances for internal nodes.}
  \item{nodesinpath }{A vector with the number of nodes along a lineage. If tipsonly==F, it will also include the number of nodes for internal nodes. These are counted from the tips to the root.}
}
\references{
Venditti, Chris, Andrew Meade, and Mark Pagel. "Detecting the node-density artifact in phylogeny reconstruction." Systematic biology 55.4 (2006): 637-643.

Hugall, Andrew F., and Michael SY Lee. "The likelihood node density effect and consequences for evolutionary studies of molecular rates." Evolution 61.10 (2007): 2293-2307.
}
\author{
David Duchene
}
\note{
Please see the references for more detailed information.
}

\seealso{
Pending.
}
\examples{
set.seed(12345)
myTree <- rtree(10)

par(mfrow = c(1, 2))
plot(myTree)
nde <- pathnode(myTree)
nde



## The function is currently defined as
function (phylo, tipsonly = T) 
{
    require(phangorn)
    di.tr <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] \%in\% phylo$edge[, 
        2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, 
        ])
    if (tipsonly == TRUE) {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == 
            root.tr, 1:length(phylo$tip.label)]
        nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo, 
            x)))
    }
    else {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == 
            root.tr, ]
	nodesinpath <- sapply(1:(length(phylo$tip.label)+phylo$Nnode), function(x) length(Ancestors(phylo, x)))
    }
    plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", 
        ylab = "Number of parent nodes", pch = 20)
    return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
  }
}
\keyword{ Node-density-effect }
