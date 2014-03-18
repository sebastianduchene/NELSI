get.lineage.time.rate <- function(taxon, tree.data.matrix){
	if(taxon %in% tree.data.matrix[, 3]){
		data.matrix <- tree.data.matrix
    	branch.times <- vector()
    	rate.time <- vector()
    	repeat{
        	parent <- data.matrix[, 2][data.matrix[, 3] == taxon]
        	time.br <- data.matrix[, 4][data.matrix[, 3] == taxon]
        	rate.br <- data.matrix[, 5][data.matrix[, 3] == taxon]
        	rate.time <- c(rate.time, rate.br)
        	branch.times <- c(branch.times, time.br)
        	taxon <- parent
        	if(!(parent %in% data.matrix[, 3])){break}
    	}
    	first.rate <- rate.time[length(rate.time)]
    	last.rate <- rate.time[1]
    	rate.time <- c(last.rate, rate.time, first.rate)
    	branch.times <- c(0, branch.times, max(tree.data.matrix[, 4]))
    	return(data.frame(branch.times, rate.time))
    }else{
    	stop("The taxon name was not found in the tree data matrix. It should be a number between 1 and the number of nodes (internal and external)")
    }
}
