#' Branch-length imbalance statistic
#'
#' Computes an imbalance statistic based on mean-centred branch lengths: the
#' sum of centred external branch lengths minus the sum of all centred branch
#' lengths, which simplifies to the negative sum of centred internal branch
#' lengths.
#'
#' @param tre A phylogenetic tree of class \code{"phylo"}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.df(tr)
#'
#' @export
get.df <-
function(tre){
	brlens <- as.numeric(scale(tre$edge.length, scale = F, center = T))
	df <- sum(brlens[which(tre$edge[,2] <= Ntip(tre))]) - sum(brlens)
	return(df)
}

