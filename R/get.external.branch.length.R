#' External (terminal) branch lengths
#'
#' Returns the lengths of all external branches of a phylogenetic tree —
#' edges whose daughter node is a tip (i.e. does not appear as a parent in
#' the edge matrix).
#'
#' @param tr A phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric vector of external branch lengths, in edge-matrix order.
#'
#' @seealso \code{\link{get.internal.branch.length}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.external.branch.length(tr)
#'
#' @export
get.external.branch.length <- function(tr){
    tr$edge.length[!tr$edge[, 2] %in% tr$edge[, 1]]
}
