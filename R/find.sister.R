#' Find tips in the sister clade
#'
#' Identifies the tip indices in the sister group of a given set of tips.
#' The sister is defined as all tips descending from the parent of the MRCA
#' of \code{clade} that are not in \code{clade}.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param clade Integer vector of tip indices whose sister group is sought.
#' @param allow.polytomy Logical. If \code{TRUE} (default), when the MRCA is
#'   a polytomy whose children are all in \code{clade} the function goes up
#'   one level to find the sister.
#'
#' @return An integer vector of tip indices constituting the sister clade.
#'
#' @seealso \code{\link{get.mrca}}, \code{\link{get.descending.nodes.branches}}
#'
#' @examples
#' library(ape)
#' tr <- read.tree(text = '((((a, b, c, d), e), f), (g, h));')
#' find.sister(tr, 1)
#'
#' @export
find.sister <- function(tr, clade, allow.polytomy = T){
tips <- 1:length(tr$tip.label)
mrca.node <- get.mrca(tr, clade)
    if(length(clade) == 1){
        all_descendants <- get.descending.nodes.branches(tr, mrca.node)$descending.nodes
        clade_tips <- all_descendants[all_descendants %in% tips]
        return(clade_tips[!clade_tips %in% clade])
    }
    if(allow.polytomy & is.polytomy(tr, mrca.node)){
        all_descendants <- get.descending.nodes.branches(tr, mrca.node)$descending.nodes
        tips_in_poly_clade <- tips[tips %in% all_descendants]
        if(all(tips_in_poly_clade %in% clade)){
            # Then go up one level
            parent_node <- tr$edge[tr$edge[, 2] == mrca.node, 1]
            nodes_from_parent_node <- get.descending.nodes.branches(tr, parent_node)$descending.nodes
            tips_from_parent_node <- tips[tips %in% nodes_from_parent_node]
            return(tips_from_parent_node[!tips_from_parent_node %in% clade])
        }else{
            return(tips_in_poly_clade[-which(tips_in_poly_clade %in% clade)])
        }
    }else{
        parent_node <- tr$edge[tr$edge[, 2] == mrca.node, 1]
        clade_nodes <- get.descending.nodes.branches(tr, mrca.node)$descending.nodes
        clade_tips <- clade_nodes[clade_nodes %in% tips]
        sister_nodes <- all_descendants[!all_descendants %in% clade_nodes]
        sister_tips <- sister_nodes[sister_nodes %in% tips]
        return(sister_tips)
    }
}


## Example:
#tr <- read.tree(text = '((((a, b, c, d), e), f), (g, h));')
#plot(tr, 'cladogram')
#nodelabels()
#tiplabels()
#find.sister(tr, 1)
#tr$tip.label
