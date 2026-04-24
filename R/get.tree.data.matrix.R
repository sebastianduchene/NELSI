#' Build the branch data matrix for a tree
#'
#' Constructs the internal matrix used by all NELSI rate-simulation functions.
#' Each row corresponds to one branch; columns record its index, parent and
#' daughter nodes, midpoint age, rate (initially \code{NA}), substitution
#' length (initially \code{NA}), and time length.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time.
#'
#' @return A numeric matrix of class \code{"tree.data.matrix"} with columns
#'   \code{branch.index}, \code{parent.node}, \code{daughter.node},
#'   \code{branch.midage}, \code{branch.rate}, \code{length.subst}, and
#'   \code{length.time}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{mid.edge.ages}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(5)
#' m <- get.tree.data.matrix(tr)
#' colnames(m)
#'
#' @export
get.tree.data.matrix <-
function(phylo){
    if (!inherits(phylo, "phylo"))
        stop("'phylo' must be an object of class \"phylo\"")
    if (!ape::is.rooted(phylo))
        stop("'phylo' must be a rooted tree")
    if (is.null(phylo$edge.length))
        stop("'phylo' must have branch lengths")
    data.matrix <- matrix(data = NA, ncol = 7, nrow = length(phylo$edge.length))
    colnames(data.matrix) <- c("branch.index", "parent.node", "daughter.node", "branch.midage", "branch.rate", "length.subst", "length.time")
    data.matrix[, 1] <- 1:length(phylo$edge.length)
    data.matrix[, 2] <- phylo$edge[ ,1]
    data.matrix[, 3] <- phylo$edge[ ,2]
    data.matrix[, 4] <- mid.edge.ages(phylo)
    data.matrix[, 7] <- phylo$edge.length
    class(data.matrix) <- "tree.data.matrix"
    return(data.matrix)
}
