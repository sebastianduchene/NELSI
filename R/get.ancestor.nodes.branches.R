#' Ancestor nodes and branches on the path to the root
#'
#' Traces the path from a target node back to the root, returning every
#' ancestor node and the edge (branch) indices traversed along the way. If
#' \code{target_node} is a vector of tip indices the function first finds their
#' MRCA and traces from there.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param target_node An integer node index, or a vector of tip indices whose
#'   MRCA is used as the starting point.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{ancestor.nodes}}{Integer vector of node indices from
#'       \code{target_node} (first) to the root (last).}
#'     \item{\code{ancestor.branches}}{Integer vector of edge row-indices in
#'       \code{tr$edge} corresponding to each step toward the root.}
#'   }
#'
#' @seealso \code{\link{get.descending.nodes.branches}}, \code{\link{get.mrca}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.ancestor.nodes.branches(tr, 3L)
#'
#' @export
get.ancestor.nodes.branches <- function(tr, target_node) {
    if (length(target_node) > 1L) {
        target_node <- castor::get_mrca_of_set(tr, target_node)
    }
    root_node <- unique(tr$edge[!(tr$edge[, 1] %in% tr$edge[, 2]), 1])
    n_clades <- length(tr$tip.label) + tr$Nnode

    # O(1) parent and edge-index lookups built once instead of scanning tr$edge each step
    parent  <- integer(n_clades)
    edge_of <- integer(n_clades)
    parent[tr$edge[, 2]]  <- tr$edge[, 1]
    edge_of[tr$edge[, 2]] <- seq_len(nrow(tr$edge))

    ancestors_nodes    <- target_node
    ancestors_branches <- integer(0)
    current <- target_node
    while (current != root_node) {
        ancestors_branches <- c(ancestors_branches, edge_of[current])
        current <- parent[current]
        ancestors_nodes <- c(ancestors_nodes, current)
    }
    list(ancestor.nodes = ancestors_nodes, ancestor.branches = ancestors_branches)
}
