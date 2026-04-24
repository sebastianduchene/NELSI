.get.label.from.pairs <- function(twolabels, useAllCharacters = FALSE) {
    if (useAllCharacters == TRUE) {
        minMax <- charMinMax(twolabels)
        k <- minMax$max
        j <- minMax$min
        newLabel <- charSum(c(charDiv2(charProd(c(k, charSum(c(k, "-1"))))), charSum(c(j, "1"))))
        return(newLabel)
    } else {
        l1 <- twolabels[1]; l2 <- twolabels[2]
        if (nchar(l1) < 14 & nchar(l2) < 14) {
            k <- max(as.numeric(l1), as.numeric(l2))
            j <- min(as.numeric(l1), as.numeric(l2))
            return(k * (k - 1) / 2 + j + 1)
        } else { return(digest::digest(sort(c(l1, l2)))) }
    }
}

.tree.labels <- function(tree, allChar = FALSE) {
    if (class(tree) != "phylo") tree <- as(tree, "phylo")
    if (is.null(tree$tip.label))
        stop("This tree has no tips")
    num.tips <- length(tree$tip.label)
    labels <- NA + 0 * (1:(nrow(tree$edge) + 1))
    names(labels)[1:num.tips] <- tree$tip.label
    names(labels)[(num.tips + 1):length(labels)] <- paste("node", 1:(length(labels) - num.tips), sep = "")
    if (allChar) {
        labels[1:num.tips] <- '1'
    } else {
        labels[1:num.tips] <- 1
    }
    NodeIDS <- (num.tips + 1):(2 * num.tips - 1)
    while (any(is.na(labels))) {
        IsReady <- NodeIDS[vapply(NodeIDS, function(x)
            !any(is.na(labels[tree$edge[which(tree$edge[, 1] == x), 2]])) & is.na(labels[x]),
            FUN.VALUE = TRUE)]
        TheseLabels <- unlist(sapply(IsReady, function(x)
            .get.label.from.pairs(labels[tree$edge[tree$edge[, 1] == x, 2]], useAllCharacters = allChar)))
        labels[IsReady] <- TheseLabels
    }
    return(labels)
}

.plot.labels <- function(tree) {
    nn <- length(tree$tip.label)
    tree$node.label <- .tree.labels(tree)[(nn + 1):(2 * nn - 1)]
    tree$tip.label  <- 1 + 0 * (1:nn)
    plot(tree, show.node.label = TRUE, edge.width = 4, cex = 1.5, edge.color = "grey")
}

.label.distance <- function(x, y) {
    uni.x   <- unique(x)
    ynotinx <- setdiff(y, x)
    dcounts <- vapply(c(uni.x, ynotinx),
                      function(k) abs(length(which(y == k)) - length(which(x == k))),
                      FUN.VALUE = 1)
    return(sum(dcounts))
}


#' Topological distance between unlabelled trees
#'
#' Computes a topological distance between two phylogenetic trees that is
#' independent of tip-label assignment. Each node is assigned a canonical
#' label derived from the labels of its descendants, and the distance is the
#' symmetric difference between the two label multisets.
#'
#' @param tree1 A rooted phylogenetic tree of class \code{"phylo"}.
#' @param tree2 A rooted phylogenetic tree of class \code{"phylo"}.
#'
#' @return A non-negative numeric scalar: the label-distance between the two
#'   trees.
#'
#' @seealso \code{\link{dist.topo.normalised}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr1 <- rcoal(6)
#' tr2 <- rcoal(6)
#' dist.unlab(tr1, tr2)
#'
#' @export
dist.unlab <- function(tree1, tree2) {
    lab1 <- .tree.labels(tree1)
    lab2 <- .tree.labels(tree2)
    return(.label.distance(lab1, lab2))
}
