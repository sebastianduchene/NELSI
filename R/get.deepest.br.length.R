#' Total length of root-descending branches
#'
#' Returns the sum of the lengths of all branches that descend directly from
#' the root node.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric scalar: the sum of root-descending branch lengths.
#'
#' @seealso \code{\link{get.oldest.branch.length}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.deepest.br.length(tr)
#'
#' @export
get.deepest.br.length <- function(tr){
  root <- tr$edge[!tr$edge[, 1] %in% tr$edge[, 2], 1][1]
  descending.branches <- tr$edge.length[tr$edge[, 1] == root]
  sum(descending.branches)
}
