require(ape)
require(phangorn)
require(geiger)

setwd("./ready_functions")
for(i in dir()) source(i)
setwd("..")

# Testing
tr1 <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 20)
s1 <- simulate.clock(tr1)

parent.daughter <- get.rate.descendant.pairs(s1[[2]])

allnode.times(s1$phylogram, tipsonly = T)

get.lineage.time.rate(4, s1$tree.data.matrix)


#

trann2trdat <- function(tree){
    require(epibase)
	tree$edge.length <- unlist(sapply(tree$annotations, function(x){ x$length }))[1:length(tree$edge.length)]
	rates <- unlist(sapply(tree$annotations, function(x){ x$rate_median }))
	if(is.ultrametric(tree) == TRUE){
		midages <- mid.edge.ages(tree)
	}else{
		midages <- mid.edge.ages(tree, max(unlist(sapply(tree$annotations, function(x){ x$height_median }))))
	}
	timelen <- tree$edge.length
	subslen <- tree$edge.length * rates
	return(data.frame(branch = rownames(as.data.frame(tree$edge)), parent = tree$edge[,1], daughter = tree$edge[,2], midage = midages, rate = rates, blensubs = subslen, blentime = timelen))
}
















get.lineage.time.rate <- function(taxon, tree.data.matrix){
	if(taxon %in% tree.data.matrix[, 3]){
		data.matrix <- tree.data.matrix
    	branch.midage <- vector()
    	rate.time <- vector()
    	repeat{
        	parent <- data.matrix[, 2][data.matrix[, 3] == taxon]
        	time.br <- data.matrix[, 4][data.matrix[, 3] == taxon]
        	rate.br <- data.matrix[, 5][data.matrix[, 3] == taxon]
        	rate.time <- c(rate.time, rate.br)
        	branch.midage <- c(branch.midage, time.br)
        	taxon <- parent
        	if(!(parent %in% data.matrix[, 3])){break}
    	}
    	first.rate <- rate.time[length(rate.time)]
    	last.rate <- rate.time[1]
    	rate.time <- c(last.rate, rate.time, first.rate)
    	branch.midage <- c(0, branch.midage, max(tree.data.matrix[, 4]))
    	return(data.frame(branch.midage, rate.time))
    }else{
    	stop("The taxon name was not found in the tree data matrix. It should be a number between 1 and the number of nodes (internal and external)")
    }
}











get.rate.descendant.pairs <- function(tree.data.matrix){
	dat <- tree.data.matrix
	parent.rate <- vector()
	daughter.rate <- vector()
	br.len <- vector()
	for(i in 1:nrow(dat)){
		daughter.temp <- dat[i, 3]
		parent.temp <- dat[i, 2]
		br.temp <- dat[i, 4]
		if(parent.temp %in% dat[, 3]){
			parent.rate <- c(parent.rate, dat[dat[, 3] == parent.temp , 5])
			daughter.rate <- c(daughter.rate, dat[i, 5])
			br.len <- c(br.len, abs(br.temp - dat[dat[, 3] == parent.temp , 5]))
		}
	}
	return(data.frame(parent.rate, daughter.rate, br.len))
}

