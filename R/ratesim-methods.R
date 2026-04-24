#' Print a rate simulation object
#'
#' Compact one-line summary of a \code{"ratesim"} object, showing tip count,
#' branch count, and the range and mean of the simulated rates.
#'
#' @param x An object of class \code{"ratesim"}.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{summary.ratesim}}, \code{\link{simulate.rate}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.clock, list(rate = 0.005, noise = 1e-6))
#' print(sim)
#'
#' @export
print.ratesim <- function(x, ...) {
    m     <- x$tree.data.matrix
    rates <- m[, "branch.rate"]
    n_tips   <- length(x$phylogram$tip.label)
    n_branches <- nrow(m)
    cat("Rate simulation (ratesim)\n")
    cat("  Tips:     ", n_tips, "\n")
    cat("  Branches: ", n_branches, "\n")
    cat("  Rates:    min =", format(min(rates),  digits = 4, scientific = TRUE),
        " | mean =",        format(mean(rates), digits = 4, scientific = TRUE),
        " | max =",         format(max(rates),  digits = 4, scientific = TRUE), "\n")
    invisible(x)
}

#' Summarise a rate simulation object
#'
#' Prints a detailed summary of a \code{"ratesim"} object including tree
#' dimensions, total branch lengths in time and substitution units, and a
#' five-number summary of the branch-rate distribution.
#'
#' @param object An object of class \code{"ratesim"}.
#' @param ... Ignored.
#'
#' @return \code{object} invisibly.
#'
#' @seealso \code{\link{print.ratesim}}, \code{\link{simulate.rate}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.uncor.lnorm,
#'                      list(mean.log = -5, sd.log = 0.5))
#' summary(sim)
#'
#' @export
summary.ratesim <- function(object, ...) {
    m <- object$tree.data.matrix
    cat("Rate simulation summary\n\n")
    cat("Tree\n")
    cat("  Tips:              ", length(object$phylogram$tip.label), "\n")
    cat("  Branches:          ", nrow(m), "\n")
    cat("  Total time length: ", round(sum(m[, "length.time"]),  6), "\n")
    cat("  Total subst length:", round(sum(m[, "length.subst"]), 6), "\n\n")
    cat("Branch rates\n")
    print(summary(m[, "branch.rate"]))
    invisible(object)
}
