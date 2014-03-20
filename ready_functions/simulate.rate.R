simulate.rate <- function(tree, FUN, ...){
	 ratesim.object <- FUN(tree, ...)
	 return(ratesim.object)
}
