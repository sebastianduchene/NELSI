\name{trann2trdat}
\alias{trann2trdat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trann2trdat(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
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
function (tree) 
{
    require(epibase)
    tree$edge.length <- unlist(sapply(tree$annotations, function(x) {
        x$length
    }))[1:length(tree$edge.length)]
    rates <- unlist(sapply(tree$annotations, function(x) {
        x$rate_median
    }))
    if (is.ultrametric(tree) == TRUE) {
        midages <- mid.edge.ages(tree)
    }
    else {
        midages <- mid.edge.ages(tree, max(unlist(sapply(tree$annotations, 
            function(x) {
                x$height_median
            }))))
    }
    timelen <- tree$edge.length
    subslen <- tree$edge.length * rates
    return(data.frame(branch = rownames(as.data.frame(tree$edge)), 
        parent = tree$edge[, 1], daughter = tree$edge[, 2], midage = midages, 
        rate = rates, blensubs = subslen, blentime = timelen))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
