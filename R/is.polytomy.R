#' Test whether a node is a polytomy
#'
#' Returns \code{TRUE} if the given node has more than two direct descendants
#' (i.e. is a polytomy / multifurcation), and \code{FALSE} if it is a
#' bifurcation or a tip.
#'
#' @param tr A phylogenetic tree of class \code{"phylo"}.
#' @param node Integer. The index of the node to test, using ape's convention:
#'   tips are \code{1..Ntips}, internal nodes are \code{Ntips+1..Ntips+Nnodes}.
#'
#' @return A logical scalar.
#'
#' @seealso \code{\link{find.monophyletic}}
#'
#' @examples
#' library(ape)
#' tr <- read.tree(text = '((a, b, c), (d, (e, f)));')
#' is.polytomy(tr, 8)   # the (a, b, c) node -- TRUE
#' is.polytomy(tr, 9)   # the root -- FALSE
#'
#' @export
is.polytomy <- function(tr, node){
  if(dim(tr$edge[tr$edge[, 1] == node, ])[1] <= 2){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

## Example:
#tr <- read.tree(text = '((a, b, c), (d, (e, f)));')
#plot(tr)
#nodelabels()
#tiplabels()

#is.polytomy(tr, 8)
#is.polytomy(tr, 9)
#is.polytomy(tr, 10)
