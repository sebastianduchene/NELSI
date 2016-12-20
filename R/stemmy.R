stemmy <- function(tre){
       if(is.ultrametric(tre)){
	stemminess <- sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))]) / sum(tre$edge.length)
       } else {
       	 tiptimes <- allnode.times(tre, tipsonly = T)
	 	  tiptimes <- as.numeric(names(tiptimes))[which(tiptimes == 0)]
		  	   stemminess <- sum(tre$edge.length[which(!tre$edge[,2] %in% tiptimes)]) / sum(tre$edge.length)
       }
       return(stemminess)
}
