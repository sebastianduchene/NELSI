get.branches.age.sorted <- function(tr){
    # Return branches sorted by the height of their terminal nodes. 
    sorted_node_heights <- sort(allnode.times(tr), dec = T)[-1]
    tr$edge.length[match(names(sorted_node_heights), tr$edge[, 2])]
}
