#' Simulate white-noise rates
#'
#' Simulates evolutionary rates along a phylogenetic tree under a white-noise
#' model. Each branch's substitution length is drawn from a lognormal
#' distribution centred on the strict-clock expectation (\code{rate *
#' branch.length}), with variance scaled by the branch's standardised
#' substitution length.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{rate}}{Numeric. The base substitution rate. Default
#'       \code{0.006}.}
#'     \item{\code{var}}{Numeric. Variance scaling parameter for the lognormal
#'       noise. Default \code{1e-6}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.clock}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.white.noise,
#'                      list(rate = 0.006, var = 1e-6))
#' plot(sim)
#'
#' @export simulate.white.noise
simulate.white.noise <-
function(tree, params = list(rate = 0.006, var = 0.000001)){
    data.matrix <- get.tree.data.matrix(tree)
    clocksubst <- tree$edge.length * params[[1]]
    clocksubstscaled <- as.numeric(scale(clocksubst)) + 1
    for(i in 1:length(clocksubst)) data.matrix[, 6][i] <- rlnorm(1, log(clocksubst[i]), abs(params[[2]] * clocksubstscaled[i]))
    branch.rates <- data.matrix[, 6] / data.matrix[, 7]
    data.matrix[, 5] <- branch.rates
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, tree.data.matrix = data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}