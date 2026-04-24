#' Root-to-tip distance versus number of ancestor nodes
#'
#' For each tip (or all nodes), computes the root-to-tip patristic distance
#' and the number of ancestral nodes on the path to the root. A scatter plot
#' of these two quantities is produced. A positive correlation indicates a
#' node-density effect, where lineages with more nodes also have longer
#' root-to-tip distances.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"}.
#' @param tipsonly Logical. If \code{TRUE} (default), only tips are included.
#'   If \code{FALSE}, all nodes (tips and internal) are included.
#'
#' @return A named list (invisibly) with two elements:
#'   \describe{
#'     \item{\code{roottotippath}}{Numeric vector of root-to-node distances.}
#'     \item{\code{nodesinpath}}{Integer vector of ancestor counts for each
#'       node.}
#'   }
#'   A scatter plot is produced as a side effect.
#'
#' @references
#' Webster, A.J., Payne, R.J.H. and Pagel, M. (2003) Molecular phylogenies
#' link rates of evolution and speciation. \emph{Science}, 301(5632), 478.
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' pdf(NULL)
#' pathnode(tr)
#' dev.off()
#'
#' @export
pathnode <- function(phylo, tipsonly = TRUE) {
    Ntips     <- length(phylo$tip.label)
    n_clades  <- Ntips + phylo$Nnode
    root_node <- phylo$edge[!(phylo$edge[, 1] %in% phylo$edge[, 2]), 1][1]

    # Build parent lookup once (O(n)) instead of scanning edges at each step
    parent <- integer(n_clades)
    parent[phylo$edge[, 2]] <- phylo$edge[, 1]

    # BFS depth (= number of ancestors) via castor's root-to-tips traversal order
    trav  <- castor::get_tree_traversal_root_to_tips(phylo, include_tips = TRUE)
    depth <- integer(n_clades)
    for (node in trav$queue) {
        if (node != root_node) depth[node] <- depth[parent[node]] + 1L
    }

    dists <- castor::get_all_distances_to_root(phylo)

    if (tipsonly) {
        idx           <- seq_len(Ntips)
        roottotippath <- dists[idx]
        names(roottotippath) <- as.character(idx)
        nodesinpath   <- depth[idx]
    } else {
        roottotippath <- dists
        names(roottotippath) <- as.character(seq_len(n_clades))
        nodesinpath   <- depth
    }

    plot(roottotippath, nodesinpath,
         xlab = "Root-to-tip path length", ylab = "Number of parent nodes", pch = 20)
    list(roottotippath = roottotippath, nodesinpath = nodesinpath)
}
