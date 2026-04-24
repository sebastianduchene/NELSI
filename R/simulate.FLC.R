#' Simulate rates under a fixed local clock model
#'
#' Simulates evolutionary rates along a phylogenetic tree under a fixed local
#' clock (FLC) model. A background rate is applied to all branches; user-
#' specified clades can override this rate with local rates, optionally
#' including or excluding their stem branches.
#'
#' @param tree A rooted chronogram of class \code{"phylo"} with branch lengths
#'   in units of time.
#' @param params A list with elements:
#'   \describe{
#'     \item{\code{clade.list}}{List of tip indices or tip labels defining the
#'       local-clock clades.}
#'     \item{\code{stem.clade.indicator}}{List of logical vectors of length 2.
#'       Each element \code{c(stem, clade)} controls whether the rate change
#'       applies to the stem branch (\code{TRUE}/\code{FALSE}) and/or to the
#'       clade branches (\code{TRUE}/\code{FALSE}).}
#'     \item{\code{background.rate}}{Numeric. Rate applied to all branches not
#'       covered by a local clock.}
#'     \item{\code{local.rates}}{List of numeric scalars, one rate per clade.}
#'   }
#'
#' @return An object of class \code{"ratesim"}; see \code{\link{simulate.rate}}.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{simulate.clock}}
#'
#' @examples
#' library(ape)
#' set.seed(1800226)
#' tr <- rcoal(10)
#' sim <- simulate.FLC(tr, params = list(
#'   clade.list            = list(c(7, 6), c(1, 2, 3)),
#'   stem.clade.indicator  = list(c(TRUE, FALSE), c(FALSE, TRUE)),
#'   background.rate       = 7e-4,
#'   local.rates           = list(28e-4, 14e-4)
#' ))
#' plot(sim)
#'
#' @export simulate.FLC
simulate.FLC <- function(tree, params = list(clade.list, stem.clade.indicator, background.rate, local.rates)){
    clade.list = params$clade.list
    stem.clade.indicator = params$stem.clade.indicator
    background.rate = params$background.rate
    local.rates = params$local.rates
        
  check.tip.labels <- function(x){
    if(is.numeric(x)){
      return(x)
    }else if(is.character(x)){
      return(match(x, tree$tip.label))
    }
  }
  
  clade.list <- lapply(clade.list, function(x) check.tip.labels(x))
  if(length(clade.list) != length(local.rates)) stop('The length of local.rates and clade.list must be thesame')
  check.monophyly <- sapply(clade.list, function(x) is.monophyletic(tree, x))
  if(any(!check.monophyly)) stop('At least one of the clades defined is not monophyletic. Please check with is.monophyletic')
  data.matrix <- get.tree.data.matrix(tree)
  data.matrix[, 'branch.rate'] <- background.rate # Note that here all branches are populated and then local clocks are populated below
  
  for(i in 1:length(clade.list)){
    # Steps below need to be done for every local clock
    mrca.node <- get.mrca(tree, clade.list[[i]])
    all.descendants <- get.descending.nodes.branches(tree, mrca.node)
    clade.branches <- all.descendants$descending.branches
    stem.branch <- data.matrix[data.matrix[, 'daughter.node'] == mrca.node, 'branch.index']
    if(stem.clade.indicator[[i]][1]){ # if it applies to the stem
      data.matrix[stem.branch, 'branch.rate'] = local.rates[[i]]
    }
    if(stem.clade.indicator[[i]][2]){ # if it applies to the clade
      data.matrix[clade.branches, 'branch.rate'] = local.rates[[i]]
    }
  }
  
    data.matrix[, 'length.subst'] <- data.matrix[, 'length.time'] * data.matrix[, 'branch.rate']
    tree$edge.length <- data.matrix[, 'length.subst']
    res <- list(tree, data.matrix)
    names(res) <- c('phylogram', 'tree.data.matrix')
    class(res) <- 'ratesim'
    return(res)
}
  

# To simulate local clock
#library(NELSI)
#set.seed(1800226)
#tree <- rcoal(10)
#par(mfrow = c(1, 1))
#plot(tree, use.edge.length = F)
#nodelabels()
#tiplabels()

#clade.list <- list(c(7, 6), c(1, 2, 3)) #"List of clades on which to set the FLCs. Can be taxon labels of indices"

#stem.clade.indicator <- list(c(T, F), c(F, T)) #"List of indicators to flag whether rate changes should be Unique to a clade and whether they should apply to the stem c(T, F), clade c(F, T), or both c(T, T)."
#background.rate <- 7e-4
#local.rates <- list(28e-4, 14e-4) # Local clock rates

#flcTest <- simulate.FLC(tree, params = list(clade.list = clade.list, 
#                                            stem.clade.indicator = stem.clade.indicator, 
#                                            background.rate = background.rate, 
#                                           local.rates = local.rates))
#plot.ratesim(flcTest, type = 's')


