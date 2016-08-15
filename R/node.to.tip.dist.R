node.to.tip.dist <-
function(phylo, tipsonly = FALSE){
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, ]
    if(tipsonly == TRUE){
    	node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, 1:length(phylo$tip.label)]
    }
    node.to.tip.distance <- abs(node.times - max(intnode.times(phylo)))
    return(node.to.tip.distance)
}
