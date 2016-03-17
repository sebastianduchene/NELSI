sum.descending.branches <- function(tr, target_node){
    tips_descendants <- match(tips(tr, target_node), tr$tip.label)
    descending_nodes <- unique(unlist(lapply(tips_descendants, function(x) nodepath(tr, from = target_node, to = x))))
    descending_nodes <- descending_nodes[descending_nodes != target_node]
    branch_lengths <- tr$edge.length[tr$edge[, 2] %in% descending_nodes]
    sum(branch_lengths)
}
