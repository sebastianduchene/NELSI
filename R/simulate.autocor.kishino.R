#' Simulate autocorrelated rates — Kishino et al. (2001) model
#'
#' Simulates evolutonary rates along a phylogenetic tree under the
#' autocorrelated lognormal model of Kishino et al. (2001). Each branch rate
#' is drawn from a lognormal distribution centred on its parent branch rate,
#' with variance proportional to branch length.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time. Branches must be ordered so that parent branches precede
#'   daughter branches (standard ape output satisfies this).
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{initial.rate}}{Numeric. Substitution rate at the root.}
#'     \item{\code{v}}{Numeric. Autocorrelation variance parameter. Larger
#'       values give weaker autocorrelation (rates diverge faster between
#'       ancestor and descendant).}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @references
#' Kishino, H., Thorne, J.L. and Bruno, W.J. (2001) Performance of a
#' divergence time estimation method under a probabilistic model of rate
#' evolution. \emph{Molecular Biology and Evolution}, 18(3), 352--361.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.autocor.thorne}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.autocor.kishino,
#'                      list(initial.rate = 0.01, v = 0.3))
#' plot(sim)
#'
#' @export simulate.autocor.kishino
simulate.autocor.kishino <- function(tree, params = list(initial.rate = 0.01, v = 0.3)) {
    require(phangorn)
    require(geiger)
    initial.rate <- params$initial.rate
    v = params$v
    data.matrix <- get.tree.data.matrix(tree)
    while(any(is.na(data.matrix[, 5])) | any(is.nan(data.matrix[, 5]))){
        data.matrix[1, 5] <- initial.rate
        for(i in 2:nrow(data.matrix)){
            parent.node <- data.matrix[i, 2]
            preceeding.parent <- data.matrix[, 2][data.matrix[ ,3] == parent.node]
            preceeding.parent.brage <- data.matrix[, 4][data.matrix[, 2] == preceeding.parent][1]
            preceeding.parent.brrate <- data.matrix[, 5][data.matrix[, 2] == preceeding.parent][1]
            if(!(is.na(preceeding.parent.brrate)) & !(is.nan(preceeding.parent.brrate)) & (parent.node %in% data.matrix[, 3])){
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), sd = v * data.matrix[i - 1, 7]^0.5))
            }else if(!(parent.node %in% data.matrix[, 3])){
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(initial.rate)), sd = sqrt(initial.rate)))
            }
        }
    }
    data.matrix[, 6]  <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
