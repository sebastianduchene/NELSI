library(NELSI)
library(ape)
source('get_tips.R')



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
     desc.tips[[i - lad_tree$Nnode - 1]] <- lad_tree$tip.label[get.tips(lad_tree, i)]
   }
   names(desc.tips) <- (1:length(desc.tips)) + length(lad_tree$tip.label)
   return(desc.tips)
}


# Next put all the trees in the correct format 
tr <- read.nexus('relaxed_uni.trees')
tr.node.mat <- get.node.mat(tr)
tr.summ.mat <- matrix(NA, nrow(tr.node.mat), 3)
tr.summ.mat[, 1] <- sapply(1:nrow(tr.node.mat), function(x) median(tr.node.mat[x, ]))
tr.summ.mat[, 2:3] <-  t(sapply(1:nrow(tr.node.mat), function(x) quantile(tr.node.mat[x, ], c(0.025, 0.975))))
colnames(tr.summ.mat) <- c('mean', 'lowerHPD', 'higherHPD')
rownames(tr.summ.mat) <- rownames(tr.node.mat)


tax.names <- tr[[1]]$tip.label

cal_raw <- xmlParse(readLines('xml_p88.txt'))
n_cals <- length(xmlChildren(xmlChildren(cal_raw)[[1]]))
cal_names <- sapply(1:n_cals, function(x) xmlAttrs(xmlChildren(cal_raw)[[1]][[x]]))

cal_list <- xmlToList(cal_raw)

node_taxa <- get.node.key(tr[[1]])

for(i in 1:length(node_taxa)){
 cat(names(node_taxa)[i], ':', file = 'relaxed_uni_names.txt', append = T)
 cat(paste(node_taxa[[i]], sep = ','), '\n', file = 'relaxed_uni_names.txt', append = T)
}


for(i in 1:length(cal_list)){
      for(k in 1:length(node_taxa)){
        if(all(node_taxa[[k]] %in% unlist(cal_list[[i]][1:(length(cal_list[[i]]) - 1)])) && all(unlist(cal_list[[i]][1:(length(cal_list[[i]]) - 1)]) %in% node_taxa[[k]]) ){
	  print(cal_list[[i]][[length(cal_list[[i]])]])
	  rownames(tr.summ.mat)[k] <- cal_list[[i]][[length(cal_list[[i]])]]
	}
      }
}

write.table(tr.summ.mat, file = 'relaxed_uni.txt')


#xmlAttrs(xmlChildren(xml_lines)[[1]][[1]][[2]])


