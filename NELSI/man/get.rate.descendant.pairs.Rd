\name{get.rate.descendant.pairs}
\alias{get.rate.descendant.pairs}
\title{
get.rate.descendant.pairs
}
\description{
The functions returns a data.frame with the rates for all pairs of parent-daughter lineage pairs. This is useful to assess the covariance of rates, which can be used as a measure of rate autcorrelation.
}
\usage{
get.rate.descendant.pairs(rate.sim.object)
}
\arguments{
  \item{rate.sim.object}{
An object of class ratesim, typically obtained with simulate.rate
}
}
\details{
None
}
\value{
A data.frame.
Each row is a daughter-parent pair. The columns are the rate of the parent branch, the rate of the daughter branch, and their difference in mid branch lengths.
}
\references{
None}
\author{
Sebastian Duchene
}
\note{
None.
}

\seealso{
plot.ratesim
ratesim
}
\examples{

set.seed(1234525)
myTree <- rcoal(50)

#Simulate rates with no autocorrelation (uncorrelated rates)

rateTreeUncor <- simulate.rate(myTree, FUN = simulate.uncor.lnorm)

uncorRates <- get.rate.descendant.pairs(rateTreeUncor)

#Simulate rates with high autocorrelation

rateTreeAutocor <- simulate.rate(myTree, FUN = simulate.autocor.kishino, params = list(initial.rate = 0.01, v = 0.01))

autocorRates <- get.rate.descendant.pairs(rateTreeAutocor)

# Plot the rates through time for all lineages to inspect the degree of autocorrelation

plot(rateTreeAutocor, col.lineages = rainbow(50))
plot(rateTreeUncor, col.lineages = rainbow(50))

# Estimate correlation of branch and daughter branch-wise rates. The correlation coefficient should be higher for the autocorrelated rates:

cor(x = autocorRates$parent.rate, y = autocorRates$daughter.rate)

cor(x = uncorRates$parent.rate, y = uncorRates$daughter.rate)


## The function is currently defined as
function (rate.sim.object) 
{
    dat <- rate.sim.object$tree.data.matrix
    parent.rate <- vector()
    daughter.rate <- vector()
    diff.br.len <- vector()
    for (i in 1:nrow(dat)) {
        daughter.temp <- dat[i, 3]
        parent.temp <- dat[i, 2]
        br.temp <- dat[i, 4]
        if (parent.temp \%in\% dat[, 3]) {
            parent.rate <- c(parent.rate, dat[dat[, 3] == parent.temp, 
                5])
            daughter.rate <- c(daughter.rate, dat[i, 5])
            diff.br.len <- c(diff.br.len, abs(br.temp - dat[dat[, 
                3] == parent.temp, 4]))
        }
    }
    return(data.frame(parent.rate, daughter.rate, diff.br.len))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ phylo }
