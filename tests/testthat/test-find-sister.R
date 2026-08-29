## Tests for find.sister()

poly_tree <- function() {
    tr <- ape::read.tree(text = "((((a, b, c, d), e), f), (g, h));")
    tr$edge.length <- rep(1, nrow(tr$edge))
    tr    # tips 1-8, nodes 9 (root) - 13
}

bif_tree <- function() {
    tr <- ape::read.tree(text = "(((a,b),(c,d)),(e,f));")
    tr$edge.length <- rep(1, nrow(tr$edge))
    tr    # tips 1-6, nodes 7 (root) - 11
}

tips_of <- function(tr, x) sort(tr$tip.label[find.sister(tr, x)])

test_that("the sister of a tip is its siblings", {
    tr <- poly_tree()
    expect_equal(tips_of(tr, 1), c("b", "c", "d"))
    expect_equal(tips_of(tr, 7), "h")
    expect_equal(tips_of(bif_tree(), 1), "b")
})

test_that("the sister of an internal node is the clade beside it, not its own tips", {
    tr <- poly_tree()
    ## node 13 is (g, h); its sister is everything else
    expect_equal(tips_of(tr, 13), c("a", "b", "c", "d", "e", "f"))
    ## node 12 is (a, b, c, d); its sister is the single tip e
    expect_equal(tips_of(tr, 12), "e")
    ## node 11 is ((a, b, c, d), e); its sister is f
    expect_equal(tips_of(tr, 11), "f")

    tr <- bif_tree()
    expect_equal(tips_of(tr, 8), c("e", "f"))   # ((a,b),(c,d))
    expect_equal(tips_of(tr, 9), c("c", "d"))   # (a,b)
})

test_that("the sister of a set of tips is found on bifurcating trees", {
    tr <- bif_tree()
    expect_equal(tips_of(tr, c(1, 2)), c("c", "d"))
    expect_equal(tips_of(tr, c(1, 2, 3, 4)), c("e", "f"))
})

test_that("polytomies are handled as documented", {
    tr <- poly_tree()
    ## part of the polytomy: the sister is the rest of it
    expect_equal(tips_of(tr, c(1, 2)), c("c", "d"))
    ## the whole polytomy: the search moves up one level
    expect_equal(tips_of(tr, c(1, 2, 3, 4)), "e")
    ## with allow.polytomy = FALSE the parent of the MRCA is always used
    expect_equal(sort(tr$tip.label[find.sister(tr, c(1, 2), allow.polytomy = FALSE)]), "e")
})

test_that("invalid input is rejected with an informative message", {
    tr <- bif_tree()
    expect_error(find.sister(tr, length(tr$tip.label) + 1L), "root has no sister")
    expect_error(find.sister(tr, 99), "must be between 1 and 11")
    expect_error(find.sister(tr, c(1, 9)), "tip indices only")
    expect_error(find.sister(tr, integer(0)), "at least one node")
    expect_error(find.sister(tr, 1.5), "node indices")
    expect_error(find.sister("not a tree", 1), "must be an object of class")
})
