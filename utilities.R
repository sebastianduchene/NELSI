##################
# UTILITIES FOR RATE SIMULATION
#####################

#####################
### FUNCTION TO DELETE SEQUENCE DATA FROM XML FILES TO OBRAIN THE MARGINAL PRIORS
#####################
get.prior.xml.file <- function(xml.file.name){
    if(xml.file.name %in% dir()){
        xml.file <- xml.file.name
        xml.name <- gsub("[.]xml", "", xml.file)
        xml.text <- readLines(xml.file)
        seq.indices <- grep("</sequence>", xml.text) - 1
        xml.text[seq.indices] <- "???"
        xml.text <- gsub(xml.name, paste0(xml.name, "_prior"), xml.text)
        new.file.name <- paste0(xml.name, "_prior", ".xml")
        writeLines(xml.text, con = new.file.name)
        return(paste0(xml.name, "_prior"))
    }else{
        stop("The file specified is not in the working directory")
    }
}
##########################################################

##################
### FUNCTION TO CONVERT THE BRANCH LENGTHS FROM A BEAST ANNOTATED TREE
##################
fix.age.annot.tree <- function(tree){
	tree$edge.length <- unlist(sapply(tree$annotations, function(x){ x$length }))[1:length(tree$edge.length)]
	return(tree)
}
##########################################################


####################
### FUNCTION TO GET THE A TREE DATA MATRIX FROM A BEAST NEXUS TREE WITH ANNOTATIONS
######################
trann2trdat <- function(tree){

	tree$edge.length <- unlist(sapply(tree$annotations, function(x){ x$length }))[1:length(tree$edge.length)]

	rates <- unlist(sapply(tree$annotations, function(x){ x$rate_median }))

	if(is.ultrametric(tree) == TRUE){

		midages <- mid.edge.ages(tree)
		
	} else {
	
		midages <- mid.edge.ages(tree, max(unlist(sapply(tree$annotations, function(x){ x$height_median }))))
	
	}
	
	timelen <- tree$edge.length
	
	subslen <- tree$edge.length * rates
	
	return(data.frame(branch = rownames(as.data.frame(tree$edge)), parent = tree$edge[,1], daughter = tree$edge[,2], midage = midages, rate = rates, blensubs = subslen, blentime = timelen))

}
##########################################################

##############
### FUNCTION TO GET THE PARENT AND DAUGHTER PAIRS OF RATES
##############
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
##############################################################

#############
## ALTERNATIVE FUNCTION TO COMPUTE THE BRANCHING TIMES OF INTERNAL AND EXTERNAL NODE REGARDLESS OF WHETHER THE TREE IS ULTRAMETRIC OR NOT
##############
allnode.times <- function(tr, tipsonly = FALSE){
    di.tr <- dist.nodes(tr)
    root.tr <- tr$edge[, 1][!(tr$edge[, 1] %in% tr$edge[, 2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, ])
    node.times <- tr.depth - di.tr[as.numeric(rownames(di.tr)) == root.tr, ]
    if(tipsonly == TRUE){
    	node.times <- tr.depth - di.tr[as.numeric(rownames(di.tr)) == root.tr, 1:length(tr$tip.label)]
    }
    return(node.times)
}
##############################################################

mid.edge.ages.2 <- function(tree){
	node.ages <- allnode.times(tree)[(length(tree$tip.label) + 1):(tree$Nnode + length(tree$tip.label))]
	tree.edges <- tree$edge
	tree.edge.lengths <- tree$edge.length
	med.times <- vector()
	for(i in 1:nrow(tree.edges)){
		med.times[i] <- (node.ages[as.numeric(names(node.ages)) == tree.edges[i,1]] - (tree.edge.lengths[i] / 2))
	}
	return(med.times)
}

####################
## FUNCTION TO GET THE MEDIAN BRANCH TIMES
################
mid.edge.ages <- function(tree){
	rootage <- max(allnode.times(tree))
	if(is.ultrametric(tree) == TRUE){

		midages <- vector()
		for(i in 1:length(tree$edge.length)){
			if(tree$edge[i, 2] > length(tree$tip.label)){
				recent.node.age <- branching.times(tree)[(tree$edge[i,2] - length(tree$tip.label))]
			halflength <- tree$edge.length[i] / 2
			midages[i] <- recent.node.age + halflength
			} else {
			midages[i] <- tree$edge.length[i] / 2
			}
		}
		return(midages)

	} else {

	nodetimes <- vector()
	extantedgelen <- max(tree$edge.length[as.vector(which(tree$edge[,1] == as.numeric(names(which(branching.times(tree) == min(branching.times(tree)))))))])
	addedval <- abs(min(branching.times(tree))) + extantedgelen
	for(i in 1:length(branching.times(tree))){
		nodetimes[i] <- (rootage / (max(branching.times(tree)) + addedval)) * (branching.times(tree) + addedval)[i]
		}

	brlen <- vector()
	for(i in 1:length(tree$edge.length)){
		brlen[i] <- (rootage / (max(branching.times(tree)) + addedval)) * tree$edge.length[i]
		}

	midages <- vector()
		for(i in 1:length(brlen)){
			if(tree$edge[i, 2] > length(tree$tip.label)){
				daughter.node.age <- nodetimes[(tree$edge[i,2] - length(tree$tip.label))]
				halflength <- brlen[i] / 2
				midages[i] <- daughter.node.age + halflength
				} else {
				parent.node.age <-  nodetimes[(tree$edge[i,1] - length(tree$tip.label))]
				midages[i] <- parent.node.age - (brlen[i] / 2)
			}
		}
		return(midages)

	}
}

##################################################################

################
## FUNCTION TO GET THE TREE DATA MATRIX
###############
get.tree.data.matrix <- function(tr){
    require(phangorn)
    require(geiger)
    data.matrix <- matrix(data = NA, ncol = 7, nrow = length(tr$edge.length))
    colnames(data.matrix) <- c("branch", "parent", "daughter", "midage", "rate", "blensubst", "blentime")
    data.matrix[, 1] <- 1:length(tr$edge.length)
    data.matrix[, 2] <- tr$edge[ ,1]
    data.matrix[, 3] <- tr$edge[ ,2]
    data.matrix[, 4] <- mid.edge.ages(tr)
    data.matrix[, 7] <- tr$edge.length
    return(data.matrix)
}
########################################################################

############
# FUNCTION TO GET THE RATE VARIATION ALONG A LINEAGE
###########
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
#######################################################

########################################
#FUNCTION TO GET THE PROPORTION OF OVERLAP BETWEEN TWO VECTORS
########################################
get.overlap <- function(dat.est, dat.sim){
    confint.est <- quantile(dat.est, c(0.025, 0.975))
    confint.sim <- quantile(dat.sim, c(0.025, 0.975))
    confint.est.width <- abs(confint.est[2] - confint.est[1])
    confint.sim.width <- abs(confint.sim[2] - confint.sim[1])
    conf.widths <- c(confint.est.width, confint.sim.width)
    bounds <- c(confint.est, confint.sim)
    if((bounds[2] > bounds[3] & bounds[2] < bounds[4])| (bounds[4] > bounds[1] & bounds[4] < bounds[2])){
        bounds.sort <- sort(bounds)
        overlap <- (bounds.sort[3] - bounds.sort[2]) / min(conf.widths)
    }else{
        overlap = 0
    }
    names(overlap) <- "overlap"
    return(overlap)
}
############################################

############################################
#FUNCTION TO GET THE MEAN RAW ERROR (PER BRANCH OR NODE)
###########################################
get.mean.error <- function(dat.est, dat.sim, nreps = 10000, type = "raw"){
    err.fun.raw <- function(x, y) return((x - y) / y)
    err.fun.abs <- function(x, y) return(abs(x - y) / y)
    if(type == "raw"){
        err.fun <- err.fun.raw
    }else if(type == "abs"){
        err.fun <- err.fun.abs
    }
    error.raw <- err.fun(dat.est, dat.sim)
    mean.error.raw <- mean(error.raw)
    mean.error.samp <- vector()
    for(i in 1:nreps){
        error.samp <- err.fun(dat.est, sample(dat.sim))
        mean.error.samp <- c(mean.error.samp, mean(error.samp))
    }
    results <- c(mean.error.raw,  1 - (sum(mean.error.raw > mean.error.samp) / nreps))
    names(results) <- c("mean.error", "random.p.val")
    return(results)
}
############################################

############################################
#FUNCTION TO GET CORRELATION BETWEEN MID-BRACH AGE AND ERROR (ABSOLUTE AND RAW)
# TAKES AS INPUT THE ESIMTATED AND SIMULATED TREE DATA MATRICES
###
###########################################
cor.error.branches <- function(dat.est, dat.sim, error.type = "raw", n.reps = 1000){
    abs.error <- function(est, sim){
        err <- abs(est - sim) / est
        return(err)
    }
    raw.error <- function(est, sim){
        err <- (est - sim) / est
        return(err)
    }
    if(error.type == "raw"){
        err.fun <- raw.error
    }else if (error.type == "abs"){
        err.fun <- abs.error
    }
    br.times <- dat.sim[, 4]
    rate.error <- err.fun(dat.est[, 5], dat.sim[, 5])
    rho <- cor(br.times, rate.error)
    rho.sims <- vector()
    for(i in 1:n.reps){
        rho.sims[i] <- cor(br.times, sample(rate.error))
    }

    results <- c(rho, sum(rho < rho.sims) / (n.reps))
    names(results) <- c("correlation", "random.p.val")
    return(results)
}
##############################################

############################################
#FUNCTION TO GET CORRELATION BETWEEN NODE AGE AND ERROR (ABSOLUTE AND RAW)
# TAKES AS INPUT THE ESIMTATED AND SIMULATED TREES
###########################################
cor.error.nodes <- function(tr.est, tr.sim, error.type = "raw", n.reps = 1000){
    error.type = "raw"
    node.ages.est <- branching.times(tr.est)
    node.ages.sim <- branching.times(tr.sim)
    node.errors.abs <- function(n.ages.est, n.ages.sim){
        error <- abs(n.ages.est - n.ages.sim) / n.ages.sim
        return(error)
    }
    node.errors.raw <- function(n.ages.est, n.ages.sim){
        error <- (n.ages.est - n.ages.sim) / n.ages.sim
        return(error)
    }
    if(error.type == "raw"){
        err.fun <- node.errors.raw
    }else if(error.type == "abs"){
        err.fun <- node.errors.abs
    }
    error.values <- err.fun(node.ages.est, node.ages.sim)
    error.cor <- cor(node.ages.est, error.values)
    error.cor.samps <- vector()
    for(i in 1:n.reps){
        error.samps.temp <- err.fun(node.ages.est, sample(node.ages.sim))
        error.cor.samps[i] <- cor(node.ages.est, error.samps.temp)
    }
    results <- c(error.cor, sum((error.cor < error.cor.samps) / (n.reps)))
    names(results) <- c("correlation", "random.p.val")
    return(results)
}
###########################################
## FUNCTION TO GET THE ROOT-TO-TIP PATH AND THE NUMBER OF NODES IN THE PATH, AND PLOT THE TWO.

pathnode <- function(phy, tipsonly = T){
	
    di.tr <- dist.nodes(phy)
    root.tr <- phy$edge[, 1][!(phy$edge[, 1] %in% phy$edge[, 2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, ])
    
    if(tipsonly == TRUE){
    	roottotippath <- di.tr[as.numeric(rownames(di.tr)) == root.tr, 1:length(phy$tip.label)]
    	nodesinpath <- sapply(1:length(phy$tip.label), function(x) length(Ancestors(phy, x)))
    } else {
		roottotippath <- di.tr[as.numeric(rownames(di.tr)) == root.tr, ]
		nodesinpath <- sapply(1:length(phy$tip.label), function(x) length(Ancestors(phy, x)))
	}
	plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", ylab = "Number of parent nodes")
	return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
}

##############################################
#FUNCTION TO CREATE HIGHLY UNBALANCED TREES. A TREE IS GENERATED WITH A COLLESS PARAMETER THAT IS WITHIN GIVEN VALUES OF STANDARD DEVIATION OF AN AVERAGE TREE OF THE SPECIFIED QUALITIES.

sim.bdtreeunbal <- function(b = 0.9, d = 0, stop = "taxa", n = 50, extinct = FALSE, sdmin = -1, sdmax = 1){
	
	cols <- vector()
	for(i in 1:100){
		cols[i] <- colless(as.treeshape(sim.bdtree(b = b, d = d, stop = stop, n = n, extinct = extinct)))
	}
	
	meancol <- mean(cols)
	sdcol <- sd(cols)
	
	minimum <- meancol + (sdmin * sdcol)
	maximum <- meancol + (sdmax * sdcol)
	
	tree <- sim.bdtree(b = b, d = d, stop = stop, n = n, extinct = extinct)
	while(minimum > colless(as.treeshape(tree)) || maximum < colless(as.treeshape(tree))){
	
	tree <- sim.bdtree(b = b, d = d, stop = stop, n = n, extinct = extinct)
	
	}
	
	return(tree)
}





