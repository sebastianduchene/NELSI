pathnode <-
function(phylo, tipsonly = T){
	 require(phangorn)
    di.tr <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, ])

    if(tipsonly == TRUE){
    	roottotippath <- di.tr[as.numeric(rownames(di.tr)) == root.tr, 1:length(phylo$tip.label)]
    	nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo, x)))
    } else {
		roottotippath <- di.tr[as.numeric(rownames(di.tr)) == root.tr, ]
		nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo, x)))
	}
	plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", ylab = "Number of parent nodes", pch = 20)
	return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
}
