\name{plot.ratesim}
\alias{plot.ratesim}
\title{
plot.ratesim
}
\description{
This function plots the rate through time for all lineages in the tree of a ratesim object in two pannels in the active graphics device. The first plot is the branch-wise rate vs the mid-age of the branch (rate vs. age). This is done for the branches that ascend from every tip in the tree. The second plot is the chronogram. The tips are coloured according to the lineages as shown in the rate vs. age plot. The width of the branches in the chronogram are proportional to their rate.
}
\usage{
plot.ratesim(rate.sim.object, col.lineages = colors(), type = "l")
}
\arguments{
  \item{rate.sim.object}{
An object of class ratesim, obtained with simulate.rate()
}
  \item{col.lineages}{
A vector with the colours for the lineages in the tree. The length should be the same as the number of tips in the tree. The default uses the first colours from the colors() function.
}
  \item{type}{
The type of plot. The default is "l", but "s" or "p"  can represent the rates more accurately because in most models the rates are averages over the branches.
}
}
\details{
none.
}
\value{
The funtion plots the rate in the active graphics device.
}
\references{
None.
}
\author{
Sebastian Duchene
}
\note{
None.
}

\seealso{
get.lineage.time.rate()
}
\examples{

set.seed(123425)
myTree <- rcoal(10)
rateTree <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot(rateTree, col.lineages = rainbow(10), type = "s")

# for a non-ultrametric tree
set.seed(123425)
myTree1 <- rtree(10)
rateTree1 <- simulate.rate(tree = myTree1, FUN = simulate.uncor.lnorm)
plot(rateTree1, col.lineages = rainbow(10), type = "s")




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
\keyword{ simulate.rate }
