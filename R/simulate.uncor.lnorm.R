#' Simulate uncorrelated lognormal rates
#'
#' Simulates evolutionary rates along a phylogenetic tree under the
#' uncorrelated lognormal (UCLN) relaxed-clock model. Each branch rate is
#' drawn independently from a lognormal distribution.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{mean.log}}{Numeric. Mean of the lognormal distribution on
#'       the log scale. Default \code{-5}.}
#'     \item{\code{sd.log}}{Numeric. Standard deviation on the log scale.
#'       Default \code{0.5}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.uncor.gamma}},
#'   \code{\link{simulate.uncor.exp}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.uncor.lnorm,
#'                      list(mean.log = -5, sd.log = 0.5))
#' plot(sim)
#'
#' @export simulate.uncor.lnorm
simulate.uncor.lnorm <-
function(tree, params = list(mean.log = -5, sd.log = 0.5)){
    mean.log <- params$mean.log
    sd.log <- params$sd.log
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rlnorm(n = length(tree$edge.length), meanlog = mean.log, sdlog = sd.log)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
