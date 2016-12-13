\name{get.ltt.summary}
\alias{get.ltt.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.ltt.summary
}
\description{
get.ltt.summary is used to obtain some summary statistics for LTT plots.
}
\usage{
get.ltt.summary(tree)
}
\arguments{
  \item{tree}{
A phylogenetic tree of class 'phylo'
}
}
\details{
This function obtains the lineages through time (LTT) for a phylogenetic tree. It uses the LTT coordinates to calculate the time of the maximum number of lineages relative to the age of the tree.
}
\value{
A vector with the ages of all nodes in the tree. The youngest tip always has age 0.
The items of the vector are numbered according to the 'phylo' object. To inspect the numbering of the tips and internal nodes see the example.
}
\references{
Some of these statisitcs have been used for ABC methods to infer epi parameters by Saulnier et al. (2016).
}
\author{
Sebastian Duchene
}
\note{
None
}
\examples{
set.seed(12345)
get_ltt_summary(rtree(100))



## The function is currently defined as

function(tree){
    if(is.ultrametric(tree)) stop("The tree is ultrametric. These statistics are only calculated for non-ultrametric \
trees")
    coords <- ltt.plot.coords(tree)
    max_location <- which.max(coords[, 2])
    relative_time_max_l <- coords[max_location, 1]/min(coords[, 1])
    reg1 <- lm(coords[1:max_location, 2] ~ coords[1:max_location, 1])
    reg2 <- lm(coords[max_location:nrow(coords), 2] ~ coords[max_location:nrow(coords), 1])
    slope1 <- reg1$coefficients[2]
    slope2 <- reg2$coefficients[2]
    ratio_slopes <- slope1/slope2
    result <- c(relative_time_max_l, slope1, slope2, ratio_slopes)
    names(result) <- c('time_max_lineages', 'slope1', 'slope2', 'ratio_slopes')
    return(result)
}

}
\keyword{ non-ultrametric tree }
