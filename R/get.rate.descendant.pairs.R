#' Parent-daughter branch rate pairs
#'
#' Extracts all ancestor-descendant branch pairs from a rate-simulation object
#' and returns their rates and the absolute difference in branch midpoint ages.
#' Useful for testing the degree of rate autocorrelation.
#'
#' @param rate.sim.object An object of class \code{"ratesim"} as returned by
#'   \code{\link{simulate.rate}}.
#'
#' @return A \code{data.frame} with columns \code{parent.rate},
#'   \code{daughter.rate}, and \code{diff.br.len} (absolute difference in
#'   branch midpoint ages).
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.autocor.kishino}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.autocor.kishino,
#'                      list(initial.rate = 0.01, v = 0.3))
#' pairs <- get.rate.descendant.pairs(sim)
#' cor(pairs$parent.rate, pairs$daughter.rate)
#'
#' @export
get.rate.descendant.pairs <-
function(rate.sim.object){
	if (!inherits(rate.sim.object, "ratesim"))
	    stop("'rate.sim.object' must be an object of class \"ratesim\"")
	dat <- rate.sim.object$tree.data.matrix
	parent.rate <- vector()
	daughter.rate <- vector()
	diff.br.len <- vector()
	for(i in 1:nrow(dat)){
		daughter.temp <- dat[i, 3]
		parent.temp <- dat[i, 2]
		br.temp <- dat[i, 4]
		if(parent.temp %in% dat[, 3]){
			parent.rate <- c(parent.rate, dat[dat[, 3] == parent.temp , 5])
			daughter.rate <- c(daughter.rate, dat[i, 5])
			diff.br.len <- c(diff.br.len, abs(br.temp - dat[dat[, 3] == parent.temp , 4]))
		}
	}
	return(data.frame(parent.rate, daughter.rate, diff.br.len))
}
