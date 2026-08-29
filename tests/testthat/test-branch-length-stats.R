## Tests for get.oldest.branch.length() and get.deepest.br.length()

## an asymmetric chronogram: the root's two children are at different ages
asym_chronogram <- function() ape::read.tree(text = "((a:1,b:1):4,(c:4,d:4):1);")

test_that("get.oldest.branch.length returns the branch of the root's oldest child", {
    tr <- asym_chronogram()
    root <- length(tr$tip.label) + 1L
    kids <- tr$edge[tr$edge[, 1] == root, 2]
    lens <- tr$edge.length[tr$edge[, 1] == root]
    ages <- max(ape::node.depth.edgelength(tr)) - ape::node.depth.edgelength(tr)

    ## it is the minimum of the root-descending branches
    expect_equal(get.oldest.branch.length(tr), min(lens))
    ## and that branch leads to the oldest of the root's children
    expect_equal(ages[kids[which.min(lens)]], max(ages[kids]))
})

test_that("get.deepest.br.length sums the root-descending branches", {
    tr <- asym_chronogram()
    root <- length(tr$tip.label) + 1L
    expect_equal(get.deepest.br.length(tr), sum(tr$edge.length[tr$edge[, 1] == root]))
})

test_that("both work on a heterochronous chronogram", {
    ## tips at different ages; the relation between branch length and age of
    ## the root's children is unchanged
    tr <- ape::read.tree(text = "((a:1,b:0.5):4,(c:4,d:3):1);")
    expect_false(ape::is.ultrametric(tr))
    root <- length(tr$tip.label) + 1L
    lens <- tr$edge.length[tr$edge[, 1] == root]
    kids <- tr$edge[tr$edge[, 1] == root, 2]
    ages <- max(ape::node.depth.edgelength(tr)) - ape::node.depth.edgelength(tr)
    expect_equal(get.oldest.branch.length(tr), min(lens))
    expect_equal(ages[kids[which.min(lens)]], max(ages[kids]))
})
