\name{get.lineage.time.rate}
\alias{get.lineage.time.rate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.lineage.time.rate(taxon, sim.rate.object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxon}{
%%     ~~Describe \code{taxon} here~~
}
  \item{sim.rate.object}{
%%     ~~Describe \code{sim.rate.object} here~~
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
function (taxon, sim.rate.object) 
{
    tree.data.matrix <- sim.rate.object[[2]]
    chrono <- sim.rate.object[[1]]
    chrono$edge.length <- tree.data.matrix[, 7]
    taxon.init <- taxon
    if (taxon \%in\% tree.data.matrix[, 3]) {
        data.matrix <- tree.data.matrix
        branch.times <- vector()
        rate.time <- vector()
        repeat {
            parent <- data.matrix[, 2][data.matrix[, 3] == taxon]
            time.br <- data.matrix[, 4][data.matrix[, 3] == taxon]
            rate.br <- data.matrix[, 5][data.matrix[, 3] == taxon]
            rate.time <- c(rate.time, rate.br)
            branch.times <- c(branch.times, time.br)
            taxon <- parent
            if (!(parent \%in\% data.matrix[, 3])) {
                break
            }
        }
        first.rate <- rate.time[length(rate.time)]
        last.rate <- rate.time[1]
        rate.time <- c(last.rate, rate.time, first.rate)
        node.times <- allnode.times(chrono)
        root.age <- max(node.times)
        branch.times <- c(node.times[taxon.init], branch.times, 
            root.age)
        return(data.frame(branch.times, rate.time))
    }
    else {
        stop("The taxon name was not found in the tree data matrix. It should be a number between 1 and the number of nodes (internal and external)")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
