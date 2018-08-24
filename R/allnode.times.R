allnode.times <- function(phylo, tipsonly = FALSE, reverse = T, keeproot = F){
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, ]
    if(reverse){
        node.times <- abs(node.times - max(node.times))
        if(keeproot){
            node.times <- node.times + phylo$root.edge
        }
    }
    if(tipsonly){
        node.times <- node.times[names(node.times) %in% 1:length(phylo$tip.label)]
    }
    return(node.times)
}
