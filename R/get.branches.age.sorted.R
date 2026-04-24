#' Branch lengths sorted by terminal-node age
#'
#' Returns branch lengths ordered from oldest to youngest terminal node. Each
#' branch is matched to its daughter node by index, so the ordering reflects
#' node age as computed by \code{\link{all.node.times}}.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time.
#'
#' @return A numeric vector of branch lengths sorted by terminal-node age
#'   (oldest first).
#'
#' @seealso \code{\link{all.node.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.branches.age.sorted(tr)
#'
#' @export
get.branches.age.sorted <- function(tr){
    # Return branches sorted by the height of their terminal nodes. 
    sorted_node_heights <- sort(all.node.times(tr), decreasing = TRUE)[-1]
    tr$edge.length[match(names(sorted_node_heights), tr$edge[, 2])]
}
