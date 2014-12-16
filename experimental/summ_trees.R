library(NELSI)
library(ape)
source('get_tips.R')

trees_sum <- lapply(1:10, function(x) rtree(20))
class(trees_sum) <- 'multiPhylo'


# First get a matrix with the ages of all internal nodes in the ladderised trees

get.node.mat <- function(tree_list){
   lad_trees <- lapply(1:length(tree_list), function(x) ladderize(tree_list[[x]]))
   time_matrix <- matrix(NA, nrow = lad_trees[[1]]$Nnode, ncol = length(lad_trees))
   for(i in 1:ncol(time_matrix)){
	time_matrix[, i] <- intnode.times(lad_trees[[i]])
   }
   rownames(time_matrix) <- (1:nrow(time_matrix)) + length(lad_trees[[1]]$tip.label)
   return(time_matrix)
}

# Get the taxa that descend from each node in the tree
get.node.key <- function(tr){
   desc.tips <- list()
   lad_tree <- ladderize(tr)
   for(i in (lad_tree$Nnode + 2):(lad_tree$Nnode+ length(lad_tree$tip.label))){
     desc.tips[[i - lad_tree$Nnode - 1]] <- get.tips(lad_tree, i)
   }
   names(desc.tips) <- (1:length(desc.tips)) + length(lad_tree$tip.label)
   return(desc.tips)
}


# Next put all the trees in the correct format 



