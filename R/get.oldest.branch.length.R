#' Length of the shortest root-descending branch
#'
#' Returns the minimum length among the branches that descend directly from
#' the root node. For a strictly bifurcating tree this is the shorter of the
#' two branches immediately below the root.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric scalar: the minimum root-descending branch length.
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
