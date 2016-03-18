get.ancestor.nodes.branches <- function(tr, target_node){
    root_node <- unique(tr$edge[!(tr$edge[, 1] %in% tr$edge[, 2]), 1])
    all_tips <- tr$edge[!(tr$edge[, 2] %in% tr$edge[, 1]), 2]
    ancestors_nodes <- target_node
    ancestors_branches <- vector()
    while(!(root_node %in% ancestors_nodes)){
        matches <- tr$edge[, 2] == ancestors_nodes[length(ancestors_nodes)]
        ancestors_nodes <- c(ancestors_nodes, tr$edge[matches, 1])
        ancestors_branches <- c(ancestors_branches, which(matches))
    }
    return(list(ancestor.nodes = ancestors_nodes, ancestor.branches = ancestors_branches))
}