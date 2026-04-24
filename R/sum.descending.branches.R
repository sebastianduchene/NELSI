#' Total branch length of a subtree
#'
#' Returns the sum of all branch lengths within the subtree rooted at
#' \code{target_node}, i.e. every branch below that node (excluding the branch
#' leading into it from above).
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param target_node Integer node index (using ape's convention:
#'   \code{1..Ntips} for tips, \code{Ntips+1..Ntips+Nnodes} for internal
#'   nodes). Passing a tip index returns 0.
#'
#' @return A single numeric value: the total branch length below
#'   \code{target_node}.
#'
#' @seealso \code{\link{get.descending.nodes.branches}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' # sum branches below root
#' sum.descending.branches(tr, length(tr$tip.label) + 1L)
#'
#' @export sum.descending.branches
sum.descending.branches <- function(tr, target_node) {
    Ntips <- length(tr$tip.label)
    if (target_node <= Ntips) return(0)
    sub <- castor::get_subtree_at_node(tr, target_node - Ntips)
    sum(sub$subtree$edge.length)
}
