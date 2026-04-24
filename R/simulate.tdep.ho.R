#' Simulate time-dependent rates — Ho et al. model
#'
#' Simulates evolutionary rates along a phylogenetic tree under a
#' time-dependent rate model. The rate at age \eqn{x} follows
#' \eqn{r(x) = s + \mu \exp(-\lambda x)}, integrating over each branch to
#' obtain the mean branch rate, with small Gaussian noise added.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{mu}}{Numeric. Amplitude of the time-dependent component.
#'       Default \code{0.035}.}
#'     \item{\code{srate}}{Numeric. Background (asymptotic) rate. Default
#'       \code{0.015}.}
#'     \item{\code{lambda}}{Numeric. Exponential decay constant. Default
#'       \code{0.1}.}
#'     \item{\code{noise}}{Numeric. Standard deviation of Gaussian noise.
#'       Default \code{0.001}.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @references
#' Ho, S.Y.W., Lanfear, R., Bromham, L., Phillips, M.J., Soubrier, J.,
#' Rodrigo, A.G. and Cooper, A. (2011) Time-dependent rates of molecular
#' evolution. \emph{Molecular Ecology}, 20(15), 3087--3101.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.clock}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.tdep.ho,
#'                      list(mu = 0.035, srate = 0.015,
#'                           lambda = 0.1, noise = 0.001))
#' plot(sim)
#'
#' @export simulate.tdep.ho
simulate.tdep.ho <-
function(tree, params = list(mu = 0.035, srate = 0.015, lambda = 0.1, noise = 0.001)){
    mu <- params$mu
    srate <- params$srate
    lambda <- params$lambda
    noise <- params$noise
    fun.rate <- function(x, m = mu, s = srate, lam = lambda){
        if(any(x >= 0)){
            return(s + (m * exp(-lam * x)))
        }else{
            stop("x is cannot be a negative number")
        }
    }

    data.matrix <- get.tree.data.matrix(tree)
    node.ages <- all.node.times(tree)
    b.times <- c(rep(0, length(tree$tip.label)), node.ages[(length(tree$tip.label) + 1):length(node.ages)])
    names(b.times) <- 1:length(b.times)

    ratetemp <- vector()
    for(i in 1:length(tree$edge.length)){
    	parentage <- b.times[as.character(data.matrix[i,2])]
    	daughterage <- b.times[as.character(data.matrix[i,3])]
    	ratetemp[i] <- integrate(fun.rate, lower = daughterage, upper = parentage)$value / data.matrix[i,7]
    }

    data.matrix[, 5] <- abs(ratetemp + rnorm(nrow(data.matrix), mean = 0, sd = noise))
    data.matrix[, 6] <- data.matrix[, 5] * data.matrix[, 7]

    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <-c("phylogram", "tree.data.matrix")
    class(res) <- "ratesim"
    return(res)
}
