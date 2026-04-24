#' Most recent common ancestor of a set of tips
#'
#' Finds the most recent common ancestor (MRCA) node of a given set of tips.
#' When a single tip is supplied the function returns its direct parent node.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param tips Integer vector of tip indices (\code{1..Ntips}).
#'
#' @return A single integer: the node index of the MRCA, following ape's
#'   convention (\code{Ntips+1..Ntips+Nnodes} for internal nodes).
#'
#' @seealso \code{\link{get.ancestor.nodes.branches}}, \code{\link{find.sister}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.mrca(tr, c(1L, 3L, 5L))
#'
#' @export
get.mrca <- function(tr, tips) {
    if (length(tips) == 1L) {
        return(tr$edge[tr$edge[, 2] == tips, 1])
    }
    castor::get_mrca_of_set(tr, descendants = tips)
}
