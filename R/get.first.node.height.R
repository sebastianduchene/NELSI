#' Age of the third-oldest internal node
#'
#' Returns the age of the third-oldest internal node (i.e. the third entry
#' when internal-node times are sorted in decreasing order). Used as a
#' summary statistic for tree shape.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"} with branch
#'   lengths in units of time.
#'
#' @return A named numeric scalar: the age of the third-oldest internal node.
#'
#' @seealso \code{\link{intnode.times}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' get.first.node.height(tr)
#'
#' @export
get.first.node.height <- function(tr){
   sort(intnode.times(tr), dec = T)[3]
}