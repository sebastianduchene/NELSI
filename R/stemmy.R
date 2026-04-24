#' Stemminess of a phylogenetic tree
#'
#' Computes stemminess as the proportion of total branch length attributable
#' to internal (non-terminal) branches. For ultrametric trees all branches
#' leading to internal nodes are counted as internal. For non-ultrametric
#' trees, branches leading to tips at time 0 (contemporary tips) are treated
#' as external.
#'
#' @param tre A phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric scalar between 0 and 1.
#'
#' @seealso \code{\link{allnode.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' stemmy(tr)
#'
#' @export
stemmy <- function(tre){
       if(is.ultrametric(tre)){
	stemminess <- sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))]) / sum(tre$edge.length)
       } else {
       	 tiptimes <- allnode.times(tre, tipsonly = T)
	 	  tiptimes <- as.numeric(names(tiptimes))[which(tiptimes == 0)]
		  	   stemminess <- sum(tre$edge.length[which(!tre$edge[,2] %in% tiptimes)]) / sum(tre$edge.length)
       }
       return(stemminess)
}
