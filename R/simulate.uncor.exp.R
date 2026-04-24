#' Simulate uncorrelated exponential rates
#'
#' Simulates evolutionary rates along a phylogenetic tree under an uncorrelated
#' relaxed-clock model where each branch rate is drawn independently from an
#' exponential distribution with the specified mean.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with element:
#'   \describe{
#'     \item{\code{mean.exp}}{Numeric. Mean of the exponential distribution.
#'       Default \code{0.001}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.uncor.lnorm}},
#'   \code{\link{simulate.uncor.gamma}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.uncor.exp, list(mean.exp = 0.001))
#' plot(sim)
#'
#' @export simulate.uncor.exp
simulate.uncor.exp <-
function(tree, params = list(mean.exp = 0.001)){
    mean.exp <- params$mean.exp
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rexp(n = length(tree$edge.length), rate = 1 / mean.exp)
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
