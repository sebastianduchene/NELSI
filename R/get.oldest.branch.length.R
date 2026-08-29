#' Length of the shortest root-descending branch
#'
#' Returns the minimum length among the branches that descend directly from
#' the root node. For a strictly bifurcating tree this is the shorter of the
#' two branches immediately below the root.
#'
#' The function is named for what that branch means on a time-scaled tree: a
#' child of the root is reached sooner along a shorter branch, so the shortest
#' root-descending branch is the one subtending the \emph{oldest} of the
#' root's children. This holds for any chronogram, ultrametric or
#' heterochronous. On a phylogram, whose branch lengths are in substitutions,
#' the returned value is simply the shortest root-descending branch and
#' carries no meaning about age, so the name should not be read literally
#' there. The tree is not checked, because a chronogram cannot be told from a
#' phylogram by inspection.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}. Give it a
#'   chronogram for the returned value to mean what the name says.
#'
#' @return A numeric scalar: the minimum root-descending branch length, that
#'   is, the length of the branch leading to the oldest child of the root.
#'
#' @seealso \code{\link{get.deepest.br.length}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.oldest.branch.length(tr)
#'
#' @export
get.oldest.branch.length <- function(tr){
    root_node <- tr$edge[!(tr$edge[, 1] %in% tr$edge[, 2]), 1][1]
    min(tr$edge.length[tr$edge[, 1] == root_node])
}
