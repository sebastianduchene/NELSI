#' Extract and plot root-to-tip distances vs sampling times
#'
#' From a rate-simulation object, extracts root-to-tip distances in
#' substitution units and in time units, produces a scatter plot, and returns
#' both as a data frame. Useful for visualising the clock signal.
#'
#' @param rate.sim.object An object of class \code{"ratesim"} as returned by
#'   \code{\link{simulate.rate}}.
#' @param tipsonly Logical. If \code{TRUE} (default), returns values for tips
#'   only.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return A \code{data.frame} with columns \code{times} (root-to-tip
#'   distances in time units) and \code{substitutions} (root-to-tip distances
#'   in substitution units). A scatter plot is produced as a side effect.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{allnode.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.uncor.lnorm,
#'                      list(mean.log = -5, sd.log = 0.5))
#' pdf(NULL)
#' get.clock.data(sim)
#' dev.off()
#'
#' @export
get.clock.data <-
function(rate.sim.object, tipsonly = T, ...){
  phylogram <- rate.sim.object$phylogram
  chrono <- rate.sim.object$phylogram
  chrono$edge.length <- rate.sim.object[[2]][, 7]
  times <- allnode.times(chrono, tipsonly)
  substitutions <- allnode.times(phylogram, tipsonly)
  plot(times, substitutions, ...)
  return(data.frame(times, substitutions))
}
