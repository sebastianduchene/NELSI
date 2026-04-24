#' Simulate strict-clock rates
#'
#' Simulates evolutionary rates along a phylogenetic tree under a strict
#' molecular clock with small Gaussian noise. All branches receive the same
#' base rate, with independent normal perturbations of standard deviation
#' \code{noise}.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{rate}}{Numeric. The clock rate. Default \code{0.006}.}
#'     \item{\code{noise}}{Numeric. Standard deviation of Gaussian noise added
#'       to each branch rate. Default \code{1e-5}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.uncor.lnorm}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.clock, list(rate = 0.005, noise = 1e-6))
#' plot(sim)
#'
#' @export simulate.clock
simulate.clock <-
function(tree, params = list(rate = 0.006, noise = 0.00001)){
    rate <- params$rate
    noise <- params$noise
    data.matrix <- get.tree.data.matrix(tree)
    branch.rates <- rep(rate, times = length(tree$edge.length))
    branch.rates <- abs(branch.rates + rnorm(length(tree$edge.length), mean = 0, sd = noise))
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
