\name{NELSI_0.2-package}
\alias{NELSI_0.2-package}
\alias{NELSI_0.2}
\docType{package}
\title{
NELSI: Nucleotide EvoLution Simulator
}
\description{
NELSI simulates rates of molecular evolution along phylogenetic trees.
}
\details{
\tabular{ll}{
Package: \tab NELSI_0.2\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2014-03-22\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
David Duchene and Sebastian Duchene
Maintainer: Sebastian Duchene <sebastian.duchene@sydney.edu.au>
}
\references{
Pending.
}
\keyword{ package }
\seealso{
}
\examples{
set.seed(1234525)

myTree <- rcoal(50)

# Simulate uncorrelated rates with default parameters:
rateTree.default <- simulate.rate(tree = myTree, FUN = simulate.uncor.lnorm)
plot(rateTree.default, col.lineages = rainbow(50))
}
