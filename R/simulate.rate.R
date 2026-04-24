#' Simulate evolutionary rates along a phylogenetic tree
#'
#' Generic wrapper that applies a rate-simulation function to a phylogenetic
#' tree. The second argument must be one of the \code{simulate.*} functions
#' provided by NELSI (or any compatible user-defined function).
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param FUN A rate-simulation function, e.g.
#'   \code{\link{simulate.uncor.lnorm}},
#'   \code{\link{simulate.autocor.kishino}},
#'   \code{\link{simulate.clock}}, etc.
#' @param ... Additional arguments passed to \code{FUN}, typically a named
#'   list \code{params = list(...)}.
#'
#' @return An object of class \code{"ratesim"} (a named list) with elements:
#'   \describe{
#'     \item{\code{phylogram}}{The input tree with branch lengths rescaled to
#'       substitution units.}
#'     \item{\code{tree.data.matrix}}{A matrix with columns
#'       \code{branch.index}, \code{parent.node}, \code{daughter.node},
#'       \code{branch.midage}, \code{branch.rate}, \code{length.subst}, and
#'       \code{length.time}.}
#'   }
#'
#' @seealso \code{\link{simulate.uncor.lnorm}},
#'   \code{\link{simulate.uncor.gamma}},
#'   \code{\link{simulate.uncor.exp}},
#'   \code{\link{simulate.autocor.kishino}},
#'   \code{\link{simulate.autocor.thorne}},
#'   \code{\link{simulate.clock}},
#'   \code{\link{simulate.dpp}},
#'   \code{\link{simulate.tdep.ho}},
#'   \code{\link{simulate.white.noise}},
#'   \code{\link{simulate.FLC}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.clock, list(rate = 0.005, noise = 1e-6))
#' class(sim)
#' names(sim)
#'
#' @export simulate.rate
simulate.rate <- function(tree, FUN, ...) {
    ratesim.object <- FUN(tree, ...)
    return(ratesim.object)
}
