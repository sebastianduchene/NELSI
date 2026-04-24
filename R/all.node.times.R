#' Ages of all nodes in a phylogenetic tree
#'
#' Returns the ages (divergence times) of all nodes — tips and internal nodes —
#' in a phylogenetic tree. Useful for analysing heterochronous trees where tip
#' ages vary.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time.
#' @param tipsonly Logical. If \code{TRUE}, returns ages of tips only. Default
#'   \code{FALSE} returns ages of all nodes.
#' @param reverse Logical. If \code{TRUE} (default), ages are expressed as time
#'   before the most recent tip (i.e. the youngest tip has age 0). If
#'   \code{FALSE}, the root has the maximum value.
#' @param keeproot Logical. If \code{TRUE} and the tree has a \code{root.edge},
#'   its length is added to all node ages. Default \code{FALSE}.
#'
#' @return A named numeric vector of node ages. Names are the integer node
#'   indices following ape convention (tips \code{1..Ntips}, internal nodes
#'   \code{Ntips+1..Ntips+Nnodes}).
#'
#' @seealso \code{\link{int.node.times}}, \code{\link{get.rtt.dist}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' all.node.times(tr)
#' all.node.times(tr, tipsonly = TRUE)
#'
#' @export all.node.times
all.node.times <- function(phylo, tipsonly = FALSE, reverse = TRUE, keeproot = FALSE) {
    if (!inherits(phylo, "phylo"))
        stop("'phylo' must be an object of class \"phylo\"")
    if (is.null(phylo$edge.length))
        stop("'phylo' must have branch lengths")
    dists <- castor::get_all_distances_to_root(phylo)
    phylo.depth <- max(dists)
    node.times <- phylo.depth - dists
    names(node.times) <- as.character(seq_along(node.times))
    if (reverse) {
        node.times <- abs(node.times - max(node.times))
        if (keeproot) node.times <- node.times + phylo$root.edge
    }
    if (tipsonly) node.times <- node.times[seq_len(length(phylo$tip.label))]
    return(node.times)
}
