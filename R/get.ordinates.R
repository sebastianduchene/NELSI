get.ordinates <- function(tr){
    nodeHeights <- allnode.times(tr, reverse = T)
    nodeHeights <- nodeHeights[!names(nodeHeights) %in% 1:length(tr$tip.label)]
    ordinates <- matrix(NA, tr$Nnode, 2)
    ordinates[, 1] <- seq(from = (tr$Nnode+length(tr$tip.label)), by = -1, length.out = nrow(ordinates))
    tips <- 1:length(tr$tip.label)
    for(i in 1:nrow(ordinates)){
        descendants <- get.descending.nodes.branches(tr, ordinates[i, 1])$descending.nodes[2:3]
        areTips <- descendants %in% tips
        if(all(areTips)){
            ordinates[i, 2] <- mean(descendants)
        }else if(sum(areTips) == 1){
            ordNode <- ordinates[ordinates[, 1] == descendants[!areTips], 2]
            ordinates[i, 2] <- mean(c(ordNode, descendants[areTips]))
        }else if(any(!areTips)){
            ordinates[i, 2] <- mean(c(ordinates[ordinates[, 1] == descendants[1], 2],
                                      ordinates[ordinates[, 1] == descendants[2], 2]))
        }
    }
    ordinates <- rbind(ordinates, cbind(length(tr$tip.label):1, length(tr$tip.label):1))
    ordinates <- ordinates[nrow(ordinates):1, ]
    ordinates <- cbind(allnode.times(tr), ordinates)
    colnames(ordinates) <- c('x.coord', 'node.index', 'y.coord')
    return(ordinates)
    #Example
    #tr <- rtree(10)
    #plot(tr, node.pos = 1, show.tip.label = F, edge.width = 4, edge.col = 'lightgrey')
    #ordinates <- get.ordinates(tr)
    #head(ordinates)
    #for(i in 1:length(tr$tip.label)){
    #    lines(rep(ordinates[i, 1], 2), c(0, ordinates[i, 3]),
    #          col = 'red', lty = 2)
    #}
    #for(i in length(tr$tip.label):nrow(ordinates)){
    #lines(rep(ordinates[i, 1], 2), c(11, ordinates[i, 3]),
    #  col = 'blue', lty = 2)
    #}
}


