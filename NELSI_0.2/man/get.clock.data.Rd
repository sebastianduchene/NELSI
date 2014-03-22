\name{get.clock.data}
\alias{get.clock.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.clock.data(rate.sim.object, tipsonly = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rate.sim.object}{
%%     ~~Describe \code{rate.sim.object} here~~
}
  \item{tipsonly}{
%%     ~~Describe \code{tipsonly} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (rate.sim.object, tipsonly = T, ...) 
{
    phylogram <- rate.sim.object$phylogram
    chrono <- rate.sim.object$phylogram
    chrono$edge.length <- rate.sim.object[[2]][, 7]
    times <- allnode.times(chrono, tipsonly)
    substitutions <- allnode.times(phylogram, tipsonly)
    plot(times, substitutions, ...)
    return(data.frame(times, substitutions))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
