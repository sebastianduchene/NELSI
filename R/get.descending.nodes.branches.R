#' Descendant nodes and branches of a subtree
#'
#' Returns all node indices and edge indices within the subtree rooted at
#' \code{target_node}, including the target node itself.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param target_node Integer node index (ape convention). Passing a tip index
#'   returns that tip only with no descendant branches.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{descending.nodes}}{Integer vector of all node indices in the
#'       subtree (tips and internal nodes), including \code{target_node}.}
#'     \item{\code{descending.branches}}{Integer vector of edge row-indices in
#'       \code{tr$edge} for every edge within the subtree.}
#'   }
#'
#' @seealso \code{\link{get.ancestor.nodes.branches}}, \code{\link{sum.descending.branches}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' root_node <- length(tr$tip.label) + 1L
#' get.descending.nodes.branches(tr, root_node)
#'
#' @export
get.descending.nodes.branches <- function(tr, target_node) {
    Ntips <- length(tr$tip.label)
    if (target_node <= Ntips) {
        return(list(descending.nodes = target_node, descending.branches = integer(0)))
    }
    sub <- castor::get_subtree_at_node(tr, target_node - Ntips)
    descending_nodes    <- c(sub$new2old_tip, sub$new2old_node + Ntips)
    descending_branches <- sub$new2old_edge
    list(descending.nodes = descending_nodes, descending.branches = descending_branches)
}
