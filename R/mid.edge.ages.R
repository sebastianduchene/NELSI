mid.edge.ages <-
function(phylo){
    require(phangorn)
	rootage <- max(allnode.times(phylo))
	if(is.ultrametric(phylo) == TRUE){

		midages <- vector()
		for(i in 1:length(phylo$edge.length)){
			if(phylo$edge[i, 2] > length(phylo$tip.label)){
				recent.node.age <- branching.times(phylo)[(phylo$edge[i,2] - length(phylo$tip.label))]
			halflength <- phylo$edge.length[i] / 2
			midages[i] <- recent.node.age + halflength
			} else {
			midages[i] <- phylo$edge.length[i] / 2
			}
		}
		return(midages)

	} else {

	nodetimes <- vector()
	extantedgelen <- max(phylo$edge.length[as.vector(which(phylo$edge[,1] == as.numeric(names(which(branching.times(phylo) == min(branching.times(phylo)))))))])
	addedval <- abs(min(branching.times(phylo))) + extantedgelen
	for(i in 1:length(branching.times(phylo))){
		nodetimes[i] <- (rootage / (max(branching.times(phylo)) + addedval)) * (branching.times(phylo) + addedval)[i]
		}

	brlen <- vector()
	for(i in 1:length(phylo$edge.length)){
		brlen[i] <- (rootage / (max(branching.times(phylo)) + addedval)) * phylo$edge.length[i]
		}

	midages <- vector()
		for(i in 1:length(brlen)){
			if(phylo$edge[i, 2] > length(phylo$tip.label)){
				daughter.node.age <- nodetimes[(phylo$edge[i,2] - length(phylo$tip.label))]
				halflength <- brlen[i] / 2
				midages[i] <- daughter.node.age + halflength
				} else {
				parent.node.age <-  nodetimes[(phylo$edge[i,1] - length(phylo$tip.label))]
				midages[i] <- parent.node.age - (brlen[i] / 2)
			}
		}
		return(round(midages, 5))

	}
}
