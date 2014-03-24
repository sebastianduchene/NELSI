\name{plot.ratesim}
\alias{plot.ratesim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.ratesim(rate.sim.object, col.lineages = colors(), type = "l")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rate.sim.object}{
%%     ~~Describe \code{rate.sim.object} here~~
}
  \item{col.lineages}{
%%     ~~Describe \code{col.lineages} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (rate.sim.object, col.lineages = colors(), type = "l") 
{
    rates.time.list <- list()
    for (i in 1:length(rate.sim.object[[1]]$tip.label)) {
        rates.time.list[[i]] <- get.lineage.time.rate(i, rate.sim.object)
    }
    ylims <- range(lapply(rates.time.list, function(y) range(y[, 
        2])))
    chrono <- rate.sim.object[[1]]
    chrono$edge.length <- rate.sim.object[[2]][, 7]
    node.ages <- allnode.times(chrono)
    xlims <- sort(range(node.ages), decreasing = T)
    par(mfrow = c(1, 2))
    plot(rates.time.list[[1]][, 1], rates.time.list[[1]][, 2], 
        ylim = ylims, xlim = xlims, ylab = "Rate", xlab = "Time", 
        type = type, lwd = 3, col = col.lineages[1])
    for (k in 2:length(rates.time.list)) {
        lines(rates.time.list[[k]][, 1], rates.time.list[[k]][, 
            2], ylim = ylims, xlim = xlims, col = col.lineages[k], 
            lwd = 3, type = type)
    }
    plot(chrono, edge.width = 1 + log(rate.sim.object[[2]][, 
        5]/min(rate.sim.object[[2]][, 5])), show.tip.label = F, 
        root.edge = T)
    tiplabels(pch = 16, col = col.lineages, cex = 1.5)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
