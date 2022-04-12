is.polytomy <- function(tr, node){
  if(dim(tr$edge[tr$edge[, 1] == node, ])[1] <= 2){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

## Example:
#tr <- read.tree(text = '((a, b, c), (d, (e, f)));')
#plot(tr)
#nodelabels()
#tiplabels()

#is.polytomy(tr, 8)
#is.polytomy(tr, 9)
#is.polytomy(tr, 10)
