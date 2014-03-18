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
