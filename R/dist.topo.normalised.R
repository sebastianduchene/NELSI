#' Normalised topological distance between two trees
#'
#' Computes the topological distance between \code{tr1} and \code{tr2} (using
#' \code{\link[ape]{dist.topo}}) and normalises it by the maximum distance
#' observed over \code{nrand} random permutations of \code{tr2}'s tip labels.
#'
#' @param tr1 A rooted phylogenetic tree of class \code{"phylo"}.
#' @param tr2 A rooted phylogenetic tree of class \code{"phylo"} with the same
#'   tip labels as \code{tr1}.
#' @param nrand Integer. Number of random permutations used to estimate the
#'   maximum distance. Default \code{100}.
#'
#' @return A numeric scalar in \code{[0, 1]}: the normalised topological
#'   distance between the two trees.
#'
#' @seealso \code{\link[ape]{dist.topo}}, \code{\link{dist.unlab}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr1 <- rcoal(10)
#' tr2 <- rcoal(10)
#' tr2$tip.label <- tr1$tip.label
#' dist.topo.normalised(tr1, tr2, nrand = 20)
#'
#' @export
dist.topo.normalised <- function(tr1, tr2, nrand=100){
    dists_sim <- vector()
    for(i in 1:nrand){
        tr_temp <- tr2
        tr_temp$tip.label <- sample(tr2$tip.label)
        dists_sim[i] <- dist.topo(tr1, tr_temp)
    }
    max_dist <- max(dists_sim)
    return(dist.topo(tr1, tr2) / max_dist)
}


