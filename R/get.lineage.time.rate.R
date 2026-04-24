#' Rate and time profile along a lineage
#'
#' Traces the path from a tip to the root and records the midpoint age and
#' rate of each branch along the lineage. The result can be used to plot how
#' the evolutionary rate changed through time on the lineage leading to a
#' given taxon.
#'
#' @param taxon Integer. Index of the tip (node index in ape convention,
#'   i.e. between 1 and the number of tips).
#' @param sim.rate.object An object of class \code{"ratesim"} as returned by
#'   \code{\link{simulate.rate}}.
#'
#' @return A \code{data.frame} with columns \code{branch.times} (age of each
#'   node along the root-to-tip path) and \code{rate.time} (branch rate at
#'   that age). The vector is padded at both ends with the root age and the
#'   tip's node time.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{plot.ratesim}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.clock, list(rate = 0.005, noise = 1e-6))
#' get.lineage.time.rate(1, sim)
#'
#' @export
get.lineage.time.rate <-
function(taxon, sim.rate.object){
        if (!inherits(sim.rate.object, "ratesim"))
            stop("'sim.rate.object' must be an object of class \"ratesim\"")
        if (!is.numeric(taxon) || length(taxon) != 1L || taxon < 1L)
            stop("'taxon' must be a single positive integer tip index")
        tree.data.matrix <- sim.rate.object[[2]]
	chrono <- sim.rate.object[[1]]
	chrono$edge.length <- tree.data.matrix[, 7]
	taxon.init <- taxon
	if(taxon %in% tree.data.matrix[, 3]){
		data.matrix <- tree.data.matrix
    	branch.times <- vector()
    	rate.time <- vector()
    	repeat{
        	parent <- data.matrix[, 2][data.matrix[, 3] == taxon]
        	time.br <- data.matrix[, 4][data.matrix[, 3] == taxon]
        	rate.br <- data.matrix[, 5][data.matrix[, 3] == taxon]
        	rate.time <- c(rate.time, rate.br)
        	branch.times <- c(branch.times, time.br)
        	taxon <- parent
        	if(!(parent %in% data.matrix[, 3])){break}
    	}
    	first.rate <- rate.time[length(rate.time)]
    	last.rate <- rate.time[1]
	rate.time <- c(last.rate, rate.time, first.rate)
	node.times <- allnode.times(chrono)
	root.age <- max(node.times)
    	branch.times <- c(node.times[taxon.init], branch.times, root.age)	
	return(data.frame(branch.times, rate.time))
    }else{
    	stop("The taxon name was not found in the tree data matrix. It should be a number between 1 and the number of nodes (internal and external)")
    }
}
