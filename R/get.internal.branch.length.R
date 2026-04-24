#' Internal branch lengths
#'
#' Returns the lengths of all internal branches of a phylogenetic tree —
#' edges whose daughter node is itself a parent (i.e. an internal node, not
#' a tip).
#'
#' @param tr A phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric vector of internal branch lengths.
#'
#' @seealso \code{\link{get.external.branch.length}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.internal.branch.length(tr)
#'
#' @export
get.internal.branch.length <- function(tr){
    internal_branches <- which(tr$edge[, 2] %in% tr$edge[, 1])
    return(tr$edge.length[internal_branches])
}
