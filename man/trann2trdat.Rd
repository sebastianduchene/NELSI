\name{trann2trdat}
\alias{trann2trdat}
\title{
trann2trdat. extract annontations from a nexus tree.
}
\description{
Convert a tree with annotations from BEAST to a data frame with the annotations for the rate, the branch lengths in substitions and the branch lengths in time. The tree should be read with funtion read.annotated.nexus from the epibase package.
}
\usage{
trann2trdat(tree)
}
\arguments{
  \item{tree}{
An object of class 'phylo' with annotations as attributes of the object. Normally the tree can be the output of BEAST or MRBayes, and it would be loaded in R with the function read.annotated.nexus from the epibase package. 
}
}
\details{
To load the tree use read.annotated.nexus("BEAST_output.tree"). This reads the annotations for each branch.
}
\value{
A data frame with the branch name, parent node, daughter node, mid age, rate, number of substitutions, and time.
}
\references{
None
}
\author{
David Duchene
}
\note{
see read.annotated.nexus from package epibase
}


\seealso{
see read.annotated.nexus from package epibase
}
\examples{
\dontrun{
myAnnotatedTree <- read.annotated.nexus("annotated.tree")
annnotationData <- trann2trdat(myAnnotatedTree)
head(annotationData)
}
## The function is currently defined as
function (tree) 
{
    require(epibase)
    tree$edge.length <- unlist(sapply(tree$annotations, function(x) {
        x$length
    }))[1:length(tree$edge.length)]
    rates <- unlist(sapply(tree$annotations, function(x) {
        x$rate_median
    }))
        midages <- mid.edge.ages(tree)
    timelen <- tree$edge.length
    subslen <- tree$edge.length * rates
    return(data.frame(branch = rownames(as.data.frame(tree$edge)), 
        parent = tree$edge[, 1], daughter = tree$edge[, 2], midage = midages, 
        rate = rates, blensubs = subslen, blentime = timelen))
  }
}
\keyword{ Phylo }
