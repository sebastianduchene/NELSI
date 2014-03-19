\name{simulate.tdep.generic}
\alias{simulate.tdep.generic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulate.tdep.generic(tree, noise = 0.001, fun = function(time) {
    0.035 * (0.015 * exp(-0.1 * time))
})
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{noise}{
%%     ~~Describe \code{noise} here~~
}
  \item{fun}{
%%     ~~Describe \code{fun} here~~
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
function (tree, noise = 0.001, fun = function(time) {
    0.035 * (0.015 * exp(-0.1 * time))
  }) 
{
    require(phangorn)
    require(geiger)
    fun.rate <- function(time) {
        if (any(time >= 0)) {
            return(fun(time))
        }
        else {
            stop("x is cannot be a negative number")
        }
    }
    data.matrix <- get.tree.data.matrix(tree)
    node.ages <- allnode.times(tree)
    b.times <- allnode.times(tree)
    b.times <- c(rep(0, length(tree$tip.label)), node.ages[(length(tree$tip.labels) + 
        1):length(node.ages)])
    names(b.times) <- 1:length(b.times)
    ratetemp <- vector()
    for (i in 1:length(tree$edge.length)) {
        parentage <- b.times[as.character(data.matrix[i, 2])]
        daughterage <- b.times[as.character(data.matrix[i, 3])]
        ratetemp[i] <- integrate(fun.rate, lower = daughterage, 
            upper = parentage)$value/data.matrix[i, 7]
    }
    data.matrix[, 5] <- abs(ratetemp + rnorm(nrow(data.matrix), 
        mean = 0, sd = noise))
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
