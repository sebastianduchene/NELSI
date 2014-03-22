\name{allnode.times}
\alias{allnode.times}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
allnode.times(phylo, tipsonly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
%%     ~~Describe \code{phylo} here~~
}
  \item{tipsonly}{
%%     ~~Describe \code{tipsonly} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
