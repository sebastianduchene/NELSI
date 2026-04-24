#' Convert an annotated tree to a branch-data matrix
#'
#' Converts a BEAST-annotated phylogenetic tree (as read by
#' \code{\link{read.annotated.nexus}} or \code{\link{read.annotated.tree}})
#' to a numeric matrix in the same format as the \code{tree.data.matrix}
#' element of a \code{"ratesim"} object.
#'
#' @param tree A phylogenetic tree of class \code{"phylo"} with a
#'   \code{$annotations} list containing \code{length} and \code{rate_median}
#'   elements for each branch.
#'
#' @return A numeric matrix with columns \code{branch}, \code{parent},
#'   \code{daughter}, \code{midage}, \code{rate}, \code{blensubs}, and
#'   \code{blentime}.
#'
#' @seealso \code{\link{read.annotated.nexus}}, \code{\link{read.annotated.tree}}
#'
#' @examples
#' \dontrun{
#' tr <- read.annotated.nexus("beast_output.trees")
#' dat <- trann2trdat(tr)
#' }
#'
#' @export
trann2trdat <-
function(tree){
	tree$edge.length <- unlist(sapply(tree$annotations, function(x){ x$length }))[1:length(tree$edge.length)]
	rates <- unlist(sapply(tree$annotations, function(x){ x$rate_median }))
	midages <- mid.edge.ages(tree)
	timelen <- tree$edge.length
	subslen <- tree$edge.length * rates
	return(cbind(branch = as.numeric(rownames(as.data.frame(tree$edge))), parent = tree$edge[,1], daughter = tree$edge[,2], midage = midages, rate = rates, blensubs = subslen, blentime = timelen))
}
