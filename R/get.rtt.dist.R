#' Root-to-tip distances for all tips
#'
#' Returns the patristic distance from the root to each tip. For ultrametric
#' trees all values are equal; for non-ultrametric (heterochronous) trees they
#' vary.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time (or substitutions per site).
#'
#' @return A named numeric vector of length \code{Ntips} with root-to-tip
#'   distances. Names are the integer tip indices (\code{1..Ntips}).
#'
#' @seealso \code{\link{node.to.tip.dist}}, \code{\link{all.node.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rtree(10)
#' get.rtt.dist(tr)
#'
#' @export
get.rtt.dist <- function(phylo) {
    dists <- castor::get_all_distances_to_root(phylo)
    Ntips <- length(phylo$tip.label)
    tip.dists <- dists[seq_len(Ntips)]
    names(tip.dists) <- as.character(seq_len(Ntips))
    return(tip.dists)
}
