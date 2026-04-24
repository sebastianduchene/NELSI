#' Root-to-tip distances for all tips
#'
#' Returns the patristic distance from the root to each tip. Equivalent to
#' \code{\link{get.rtt.dist}}.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time (or substitutions per site).
#'
#' @return A named numeric vector of length \code{Ntips} with root-to-tip
#'   distances. Names are the integer tip indices (\code{1..Ntips}).
#'
#' @seealso \code{\link{get.rtt.dist}}, \code{\link{all.node.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rtree(10)
#' node.to.tip.dist(tr)
#'
#' @export
node.to.tip.dist <- function(phylo) {
    dists <- castor::get_all_distances_to_root(phylo)
    Ntips <- length(phylo$tip.label)
    tip.dists <- dists[seq_len(Ntips)]
    names(tip.dists) <- as.character(seq_len(Ntips))
    return(tip.dists)
}
