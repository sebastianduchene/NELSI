#' Midpoint ages of all branches
#'
#' Computes the age at the midpoint of each branch in a phylogenetic tree.
#' For ultrametric trees the midpoint is the daughter-node age plus half the
#' branch length. For non-ultrametric trees the ages are first rescaled so
#' that the root age equals the maximum \code{allnode.times} value, then the
#' same formula is applied.
#'
#' @param phylo A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time.
#'
#' @return A numeric vector of midpoint ages, one per branch, in the same
#'   order as \code{phylo$edge}.
#'
#' @seealso \code{\link{get.tree.data.matrix}}, \code{\link{allnode.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(5)
#' mid.edge.ages(tr)
#'
#' @export
mid.edge.ages <-
function(phylo){
	rootage <- max(allnode.times(phylo))
	if(is.ultrametric(phylo) == TRUE){

		midages <- vector()
		for(i in 1:length(phylo$edge.length)){
			if(phylo$edge[i, 2] > length(phylo$tip.label)){
				recent.node.age <- branching.times(phylo)[(phylo$edge[i,2] - length(phylo$tip.label))]
			halflength <- phylo$edge.length[i] / 2
			midages[i] <- recent.node.age + halflength
			} else {
			midages[i] <- phylo$edge.length[i] / 2
			}
		}
		return(midages)

	} else {

	nodetimes <- vector()
	extantedgelen <- max(phylo$edge.length[as.vector(which(phylo$edge[,1] == as.numeric(names(which(branching.times(phylo) == min(branching.times(phylo)))))))])
	addedval <- abs(min(branching.times(phylo))) + extantedgelen
	for(i in 1:length(branching.times(phylo))){
		nodetimes[i] <- (rootage / (max(branching.times(phylo)) + addedval)) * (branching.times(phylo) + addedval)[i]
		}

	brlen <- vector()
	for(i in 1:length(phylo$edge.length)){
		brlen[i] <- (rootage / (max(branching.times(phylo)) + addedval)) * phylo$edge.length[i]
		}

	midages <- vector()
		for(i in 1:length(brlen)){
			if(phylo$edge[i, 2] > length(phylo$tip.label)){
				daughter.node.age <- nodetimes[(phylo$edge[i,2] - length(phylo$tip.label))]
				halflength <- brlen[i] / 2
				midages[i] <- daughter.node.age + halflength
				} else {
				parent.node.age <-  nodetimes[(phylo$edge[i,1] - length(phylo$tip.label))]
				midages[i] <- parent.node.age - (brlen[i] / 2)
			}
		}
		return(round(midages, 5))

	}
}
