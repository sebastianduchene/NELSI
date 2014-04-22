
read.beast.tree <- function(file){
	tree <- read.annotated.nexus(file)
	data.matrix <- trann2trdat(tree)		
	res <- list(tree, data.matrix)
	names(res) <- c("chronogram", "tree.data.matrix")
	class(res) <- "beast.tree"
	return(res)
}

tre <- read.beast.tree(file)
plot(tre)
