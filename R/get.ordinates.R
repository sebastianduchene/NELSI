#' Node plotting coordinates for a phylogenetic tree
#'
#' Computes x and y plotting coordinates for every node (tips and internal
#' nodes) of a phylogenetic tree. x coordinates are node ages (from
#' \code{\link{all.node.times}}); y coordinates are based on the average
#' ordinate of direct descendants, assigned from leaves upward. Used
#' internally by \code{\link{plot.tree.lines}}.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric matrix with columns \code{x.coord} (node age),
#'   \code{node.index} (ape node index), and \code{y.coord} (plotting
#'   ordinate).
#'
#' @seealso \code{\link{plot.tree.lines}}, \code{\link{all.node.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' head(get.ordinates(tr))
#'
#' @export
get.ordinates <- function(tr) {
    Ntips <- length(tr$tip.label)

    # Post-order traversal: process internal nodes leaves-to-root so every
    # node's y.coord is computed before its parent needs it.
    trav       <- castor::get_tree_traversal_root_to_tips(tr, include_tips = FALSE)
    post_order <- rev(trav$queue)   # queue already uses ape node indices

    y_coord <- numeric(Ntips + tr$Nnode)
    y_coord[seq_len(Ntips)] <- seq_len(Ntips)   # tips: y = tip index (1..Ntips)

    for (node in post_order) {
        children       <- tr$edge[tr$edge[, 1] == node, 2]
        y_coord[node]  <- mean(y_coord[children])
    }

    x_coord  <- all.node.times(tr)
    node_idx <- seq_len(Ntips + tr$Nnode)
    ordinates <- cbind(x.coord = x_coord, node.index = node_idx, y.coord = y_coord)
    return(ordinates)
}
