allnode.times <-
function(phylo, tipsonly = FALSE, reverse = T, keeproot = F){
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, ]
#####
    if((!is.null(phylo$root.edge)) & keeproot){
        node.times <- (max(node.times)+tr$root.edge) - node.times
    }
#####
    if(tipsonly){
        node.times <- node.times[names(node.times) %in% 1:length(phylo$tip.label)]
    }
#    if(tipsonly == TRUE){
#    	node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, 1:length(phylo$tip.label)]
#    }
    if(reverse){
        node.times <- abs(node.times - max(node.times))
    }
    return(node.times)
}
