\name{mid.edge.ages}
\alias{mid.edge.ages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mid.edge.ages(phylo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
%%     ~~Describe \code{phylo} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
