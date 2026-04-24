#' Ages of internal nodes in a phylogenetic tree
#'
#' Returns the ages (divergence times) of internal nodes only, expressed as
#' time before the most recent tip. Equivalent to the forward-time distances
#' from root used internally by other NELSI functions.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time.
#'
#' @return A named numeric vector of length \code{Nnodes}. Names are the
#'   integer internal-node indices (\code{Ntips+1..Ntips+Nnodes}).
#'
#' @seealso \code{\link{all.node.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' int.node.times(tr)
#'
#' @export
int.node.times <- function(phylo) {
    dists <- castor::get_all_distances_to_root(phylo)
    phylo.depth <- max(dists)
    int.idx <- seq(length(phylo$tip.label) + 1L, length(dists))
    node.times <- phylo.depth - dists[int.idx]
    names(node.times) <- as.character(int.idx)
    return(node.times)
}
