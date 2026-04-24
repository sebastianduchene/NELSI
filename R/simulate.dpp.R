#' Simulate rates under a Dirichlet process prior
#'
#' Simulates evolutionary rates along a phylogenetic tree under a Dirichlet
#' process prior (DPP) model. Branches are grouped into rate categories by
#' sampling from a DPP, and each category is assigned a rate drawn from a
#' gamma distribution.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{alpha}}{Numeric. Concentration parameter of the Dirichlet
#'       process. Larger values produce more rate categories. Default \code{1}.}
#'     \item{\code{shape}}{Numeric. Shape parameter of the gamma distribution
#'       for rate categories. Default \code{3.98}.}
#'     \item{\code{rate}}{Numeric. Rate parameter of the gamma distribution.
#'       Default \code{516.53}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.uncor.gamma}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.dpp,
#'                      list(alpha = 1, shape = 3.98, rate = 516.53))
#' plot(sim)
#'
#' @export simulate.dpp
simulate.dpp <-
function(tree, params = list(alpha = 1, shape = 3.98, rate = 516.53)){
    rdirichlet <- function(n, alpha) {
        x <- matrix(rgamma(n * length(alpha), shape = alpha, rate = 1),
                    nrow = n, byrow = TRUE)
        x / rowSums(x)
    }
    shape.gamma <- params$shape
    rate.gamma <- params$rate
    alpha <- params$alpha
    data.matrix <- get.tree.data.matrix(tree)
    nbranches <- nrow(tree$edge)
    cats <- sample(1:nbranches, 1, prob = rdirichlet(1, rep(alpha, nbranches)))
    branch.cats <- sample(1:cats, nbranches, replace = T, prob = rdirichlet(1, rep(alpha, cats)))
    branch.rates <- rgamma(n = cats, shape = shape.gamma, rate = rate.gamma)
    names(branch.rates) <- 1:cats
    branch.rates <- branch.rates[branch.cats]
    data.matrix[, 5] <- branch.rates
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
