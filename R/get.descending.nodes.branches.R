get.descending.nodes.branches <- function(tr, target_node){
    root_node <- unique(tr$edge[!(tr$edge[, 1] %in% tr$edge[, 2]), 1])
    all_tips <- tr$edge[!(tr$edge[, 2] %in% tr$edge[, 1]), 2]
    descendant_nodes <- target_node 
    descendant_branches <- vector()
    nodes_temp <- target_node
    while(!all(nodes_temp %in% all_tips)){
        matches <- tr$edge[, 1] %in% nodes_temp
        nodes_temp <- tr$edge[matches, 2]
        descendant_nodes <- c(descendant_nodes, nodes_temp)
        descendant_branches <- c(descendant_branches, which(matches))
    }
    return(list(descending.nodes = descendant_nodes, descending.branches = descendant_branches))
}