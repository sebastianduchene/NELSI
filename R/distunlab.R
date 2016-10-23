require(ape)
require(digest)

getLabelFromPairs <-function( twolabels,useAllCharacters=FALSE ) {
    if (useAllCharacters==TRUE) {
        minMax=charMinMax(twolabels)
	k = minMax$max
	j = minMax$min
	newLabel=charSum(c(charDiv2(charProd(c(k, charSum(c(k,"-1"))))),charSum(c(j,"1")))) #this is k*(k+1)/2 + j + 1 (include the 1 if 0 is a tip)
        return(newLabel)
    } else {
        l1=twolabels[1]; l2=twolabels[2];
        if (nchar(l1) < 14 & nchar(l2) < 14) {
            k = max(as.numeric(l1),as.numeric(l2))
            j = min(as.numeric(l1),as.numeric(l2))
            return( k*(k-1)/2 + j + 1) # NOTE if 1 is a tip and there are no 0s allowed (full binary tree) the correct expression is 1/2 k (k-1) + j + 1.
        } else { return(digest(sort(c(l1,l2)))) }
    }
}

treelabels <- function(tree,allChar=FALSE) {
    if (class(tree) != "phylo") tree = as(tree,"phylo")
     if (is.null(tree$tip.label))
         stop("This tree has no tips")
    num.tips=length(tree$tip.label)
    labels=NA + 0*(1:(nrow(tree$edge)+1))
    names(labels)[1:num.tips]=tree$tip.label;
    names(labels)[(num.tips+1): length(labels)]=paste("node",1:(length(labels)-num.tips),sep="")
    if(allChar){
        labels[1:num.tips] <- '1'
    }else{
        labels[1:num.tips] <- 1
    }
    NodeIDS= (num.tips + 1) : (2*num.tips -1)
    while (any(is.na(labels))) {
        IsReady = NodeIDS[ vapply(NodeIDS,function(x) !any(is.na(labels[tree$edge[which(tree$edge[,1]==x),2]])) & is.na(labels[x])  ,FUN.VALUE=TRUE) ]
        TheseLabels = unlist(sapply(IsReady, function(x) getLabelFromPairs(labels[tree$edge[tree$edge[,1]==x,2]],useAllCharacters=allChar)))
        labels[IsReady]=TheseLabels
    }
    return(labels)
}

plotlabels <- function(tree) {
    nn=length(tree$tip.label)
    tree$node.label=treelabels(tree)[(nn+1):(2*nn-1)]
    tree$tip.label=1+0*(1:nn)
    plot(tree,show.node.label=TRUE, edge.width=4,cex=1.5,edge.color="grey")
}


labeldistance <- function(x,y) {
	# for each unique element of x, how many times does it come up in x, and how many in y?
    uni.x=unique(x)
    ynotinx=setdiff(y,x) # things in y not in x
    dcounts=vapply(c(uni.x,ynotinx), function(k) abs(length(which(y==k))-length(which(x==k))),FUN.VALUE=1)
    return(sum(dcounts))
	# for each unique element of y NOT already counted in x, how many times does it happen in y?
}	# distance is the sum of all of those numbers.



distunlab<-function(tree1,tree2) {
    lab1=treelabels(tree1)
    lab2=treelabels(tree2)
	return(labeldistance(lab1,lab2))
}
