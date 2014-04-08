simulate.autocor.thorne <-
function(tree, params = list(initial.rate = 0.01, v = 0.3)){
    require(phangorn)
    require(geiger)
    initial.rate <- params$initial.rate
    v = params$v
    data.matrix <- get.tree.data.matrix(tree)
    while(any(is.na(data.matrix[, 5])) | any(is.nan(data.matrix[, 5]))){
        data.matrix[1, 5] <- initial.rate
        for(i in 2:nrow(data.matrix)){
            parent.node <- data.matrix[i, 2]
            preceeding.parent <- data.matrix[, 2][data.matrix[ ,3] == parent.node]
            preceeding.parent.brage <- data.matrix[, 4][data.matrix[, 2] == preceeding.parent][1]
            preceeding.parent.brrate <- data.matrix[, 5][data.matrix[, 2] == preceeding.parent][1]
            if(!(is.na(preceeding.parent.brrate)) & !(is.nan(preceeding.parent.brrate)) & (parent.node %in% data.matrix[, 3])){
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), sd = v * abs(data.matrix[i, 4] - preceeding.parent.brage)^0.5))
            }else if(!(parent.node %in% data.matrix[, 3])){
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(initial.rate)), sd = sqrt(initial.rate)))
            }
        }
    }
    data.matrix[, 6]  <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    class(res)  <- "ratesim"
    return(res)
}
