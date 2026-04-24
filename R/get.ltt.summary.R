#' Lineage-through-time summary statistics
#'
#' Summarises the lineage-through-time (LTT) curve of a non-ultrametric tree
#' using two linear regressions: one before and one after the peak lineage
#' count. Only valid for non-ultrametric (heterochronous) trees.
#'
#' @param tree A non-ultrametric phylogenetic tree of class \code{"phylo"}.
#'   An error is raised if the tree is ultrametric.
#'
#' @return A named numeric vector of length 4:
#'   \describe{
#'     \item{\code{time_max_lineages}}{Relative time (fraction of total tree
#'       age) at which the maximum number of lineages occurs.}
#'     \item{\code{slope1}}{Slope of the linear regression before the
#'       diversity peak (increasing phase).}
#'     \item{\code{slope2}}{Slope of the linear regression after the diversity
#'       peak (decreasing phase).}
#'     \item{\code{ratio_slopes}}{Ratio \code{slope1 / slope2}.}
#'   }
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rtree(20)
#' get.ltt.summary(tr)
#'
#' @export
get.ltt.summary <- function(tree) {
    if (is.ultrametric(tree)) stop("The tree is ultrametric. These statistics are only calculated for non-ultrametric trees")
    coords <- ltt.plot.coords(tree)
    max_location <- which.max(coords[, 2])
    relative_time_max_l <- coords[max_location, 1] / min(coords[, 1])
    reg1 <- lm(coords[1:max_location, 2] ~ coords[1:max_location, 1])
    reg2 <- lm(coords[max_location:nrow(coords), 2] ~ coords[max_location:nrow(coords), 1])
    slope1 <- reg1$coefficients[2]
    slope2 <- reg2$coefficients[2]
    ratio_slopes <- slope1 / slope2
    result <- c(relative_time_max_l, slope1, slope2, ratio_slopes)
    names(result) <- c('time_max_lineages', 'slope1', 'slope2', 'ratio_slopes')
    return(result)
}
