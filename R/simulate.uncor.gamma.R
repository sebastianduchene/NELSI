#' Simulate uncorrelated gamma-distributed rates
#'
#' Simulates evolutionary rates along a phylogenetic tree under an uncorrelated
#' relaxed-clock model where each branch rate is drawn independently from a
#' gamma distribution. If \code{rate} is \code{NULL} it is inferred as
#' \code{shape / mean}.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{mean}}{Numeric. Mean of the gamma distribution. Default
#'       \code{0.007}.}
#'     \item{\code{shape}}{Numeric. Shape parameter. Default \code{3.98}.}
#'     \item{\code{rate}}{Numeric or \code{NULL}. Rate parameter. If
#'       \code{NULL} (default), computed as \code{shape / mean}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.uncor.lnorm}},
#'   \code{\link{simulate.uncor.exp}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.uncor.gamma,
#'                      list(mean = 0.007, shape = 3.98))
#' plot(sim)
#'
#' @export simulate.uncor.gamma
simulate.uncor.gamma <-
function(tree, params = list(mean = 0.007, shape = 3.98, rate = NULL)){
    if(is.null(params$rate)) params$rate <- params$shape/params$mean
    rate.gamma <- params$rate
    shape.gamma <- params$shape
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rgamma(n = length(tree$edge.length), shape = shape.gamma, rate = rate.gamma)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
