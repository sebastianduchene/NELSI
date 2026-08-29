#' Find tips in the sister clade
#'
#' Identifies the tips in the sister group of a node, or of a set of tips.
#' The sister group of a node is made up of the tips that descend from its
#' parent but not from the node itself.
#'
#' \code{clade} may be a single node, either a tip or an internal node, or a
#' set of two or more tips, in which case the sister group of their most
#' recent common ancestor is returned. When a set of tips does not account for
#' every tip below their MRCA and that MRCA is a polytomy, the remaining tips
#' of the polytomy are returned instead (see \code{allow.polytomy}).
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param clade Either a single node index (a tip, \code{1..Ntip}, or an
#'   internal node, \code{Ntip+1..Ntip+Nnode}), or an integer vector of two or
#'   more tip indices.
#' @param allow.polytomy Logical. If \code{TRUE} (default), and \code{clade}
#'   is a set of tips whose MRCA is a polytomy, the sister is taken to be the
#'   other tips of that polytomy; when \code{clade} accounts for all of them,
#'   the search moves up one level. If \code{FALSE}, the sister is always
#'   taken from the parent of the MRCA.
#'
#' @return An integer vector of tip indices, in increasing order, constituting
#'   the sister clade.
#'
#' @seealso \code{\link{get.mrca}}, \code{\link{get.descending.nodes.branches}}
#'
#' @examples
#' library(ape)
#' tr <- read.tree(text = '((((a, b, c, d), e), f), (g, h));')
#' find.sister(tr, 1)     # sister of tip a
#' find.sister(tr, 13)    # sister of the internal node (g, h)
#' find.sister(tr, c(1, 2, 3, 4))
#'
#' @export
find.sister <- function(tr, clade, allow.polytomy = TRUE){
    if (!inherits(tr, "phylo"))
        stop("'tr' must be an object of class \"phylo\"")
    n.tip <- length(tr$tip.label)
    n.total <- n.tip + tr$Nnode
    root <- n.tip + 1L

    if (!length(clade))
        stop("'clade' must contain at least one node")
    if (!is.numeric(clade) || anyNA(clade) || any(clade %% 1 != 0))
        stop("'clade' must be given as node indices")
    clade <- as.integer(clade)
    if (any(clade < 1L) || any(clade > n.total)) {
        stop("'clade' must be between 1 and ", n.total,
             " (tips are 1..", n.tip, ", internal nodes ", root, "..", n.total, ")")
    }

    ## tips descending from a node, the node itself included when it is a tip
    descending.tips <- function(node) {
        d <- get.descending.nodes.branches(tr, node)$descending.nodes
        sort(d[d <= n.tip])
    }

    ## the node whose sister group is wanted
    if (length(clade) == 1L) {
        node <- clade
    } else {
        if (any(clade > n.tip)) {
            stop("when 'clade' has more than one element it must contain tip ",
                 "indices only (1..", n.tip, "); node ",
                 paste(clade[clade > n.tip], collapse = ", "), " is internal")
        }
        node <- get.mrca(tr, clade)
    }
    if (node == root)
        stop("the root has no sister clade")

    node.tips <- descending.tips(node)

    ## a set of tips that does not account for the whole polytomy below their
    ## MRCA: the sister is the rest of that polytomy
    if (length(clade) > 1L && allow.polytomy && is.polytomy(tr, node) &&
        !all(node.tips %in% clade)) {
        return(setdiff(node.tips, clade))
    }

    parent <- tr$edge[tr$edge[, 2] == node, 1]
    setdiff(descending.tips(parent), node.tips)
}
