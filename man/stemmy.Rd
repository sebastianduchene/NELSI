\name{stemmy}
\alias{stemmy}
\title{
stemmy
}
\description{
stemmy calculates the stemminess of a phylogenetic tree. 
}
\usage{
stemmy(tre)
}
\arguments{
  \item{tre}{
An object of class phylo.
}
}
\details{
None.
}
\value{
numeric.
}

\author{
David Duchene
}



## The function is currently defined as
function(tre)
{
       if(is.ultrametric(tre)){
	stemminess <- sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))]) / sum(tre$edge.length)
       } else {
       	 tiptimes <- allnode.times(tre, tipsonly = T)
	 	  tiptimes <- as.numeric(names(tiptimes))[which(tiptimes == 0)]
		  	   stemminess <- sum(tre$edge.length[which(!tre$edge[,2] %in% tiptimes)]) / sum(tre$edge.length)
       }
       return(stemminess)
}
\keyword{ molecular-clock }
