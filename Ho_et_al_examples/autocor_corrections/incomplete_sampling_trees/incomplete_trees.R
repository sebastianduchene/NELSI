library(phangorn)
library(phangorn)
library(geiger)
library(NELSI)


root_tips <- function(tr){
  n_tips <- length(tr$tip.label)
  tips_index <- 1:n_tips
  root_index <- n_tips + 1
  root_tips <- tips_index[tips_index %in% Children(tr, root_index)]
  if(length(root_tips) == 0){
    root_tips <- NULL
  }
  return(root_tips)
}

rem_tips <- function(tr){
  tips_index <- tr$tip.label[!(1:length(tr$tip.label) %in% root_tips(tr))]
}


prune_tree_random <- function(tr, n){
  tr_pruned <- tr
  if(n > (length(tr$tip.label) - 2)) stop("The number of tips should not be larger than the number of tips - 2")
  all_rem_tips <- vector()
  for(i in 1:n){
    removable_tips <- rem_tips(tr_pruned)
    tips_removed <- sample(removable_tips, 1)
    tr_pruned <- ape::drop.tip(tr_pruned, tips_removed)
    all_rem_tips <- c(all_rem_tips, tips_removed)
  }
  return(all_rem_tips)
}


#################################
# Testing the functions
#################################
#tr_test <- sim.bdtree(b = 1, d = 0, stop = 'taxa', n = 10, extinct = FALSE)

#tr_test$edge.length <- tr_test$edge.length * (10 / max(allnode.times(tr_test)))

# To prune the tree: Write a function to get the tips that directly descend from the root. These are excluded from the sampling each time.


#plot(tr_test, show.tip.label = F)
#tiplabels()
#nodelabels()

#keep_tips_t1 <- root_tips(tr_test)
#rem_tips_t1 <- rem_tips(tr_test)

#print(paste("the root tip is:", keep_tips_t1))
#print(paste("The removable tips are:", paste(rem_tips_t1, collapse = ", ")))

#print("Example pruning 5 tips from the tree")
#print(tr_test)
#print(prune_tree_random(tr_test, 5))

##################################
# Simulate ten trees with 100 taxa and remove 50
##################################

#for(i in 1:10){
#      tr_temp_comp <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = 100, extinct = F)
#      tr_temp_comp$edge.length <- tr_temp_comp$edge.length * ( 10 / max(allnode.times(tr_temp_comp)))
#      tr_temp_prune <- prune_tree_random(tr_temp_comp, 50)
#      write.tree(tr_temp_comp, file = paste0("comp_yule_", i, ".tre"))
#      write.tree(tr_temp_prune, file = paste0("pruned_yule_", i, ".tre"))
#}
