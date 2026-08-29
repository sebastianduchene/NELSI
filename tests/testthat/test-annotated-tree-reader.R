## Tests for read.annotated.nexus()

example_tree <- function() {
    for (p in c("example_data/hiv_A_env.tree",
                "../../example_data/hiv_A_env.tree",
                testthat::test_path("..", "..", "example_data", "hiv_A_env.tree"))) {
        if (file.exists(p)) return(p)
    }
    testthat::skip("example_data/hiv_A_env.tree not found")
}

## A small NEXUS file holding several annotated trees, written on the fly.
multi_tree_file <- function(n = 3) {
    f <- tempfile(fileext = ".trees")
    trees <- vapply(seq_len(n), function(i) sprintf(
        "\ttree STATE_%d = ((1[&rate=0.1]:%f,2[&rate=0.2]:%f)[&rate=0.3]:0.5,3[&rate=0.4]:1.0);",
        i, 0.5 + i / 100, 0.5 + i / 100), character(1))
    writeLines(c("#NEXUS", "", "Begin trees;", "\tTranslate",
                 "\t\t1 taxon_a,", "\t\t2 taxon_b,", "\t\t3 taxon_c", "\t\t;",
                 trees, "End;"), f)
    f
}

test_that("a single annotated tree is read with the correct tip labels", {
    f <- example_tree()
    tr <- read.annotated.nexus(f)
    ref <- ape::read.nexus(f)

    expect_s3_class(tr, "phylo")
    ## the tips of this file are NOT in TRANSLATE-block order in the Newick
    ## string, which used to produce permuted tip labels
    expect_true(isTRUE(ape::all.equal.phylo(tr, ref, use.edge.length = TRUE)))
    expect_setequal(tr$tip.label, ref$tip.label)
})

test_that("annotations are named by node, following the row order of $edge", {
    tr <- read.annotated.nexus(example_tree())
    root <- length(tr$tip.label) + 1L

    ## one annotation per edge, plus the root
    expect_equal(length(tr$annotations), nrow(tr$edge) + 1L)
    expect_identical(names(tr$annotations),
                     as.character(c(tr$edge[, 2], root)))
    expect_false(anyDuplicated(names(tr$annotations)) > 0)
    ## every node of the tree is covered exactly once
    expect_setequal(as.integer(names(tr$annotations)),
                    seq_len(length(tr$tip.label) + tr$Nnode))
    ## and can be reached by node number
    expect_true(is.list(tr$annotations[[as.character(root)]]))
})

## A single tree whose tips appear in the Newick string in an order that is NOT
## the TRANSLATE-block order (3, 1, 4, 2), each node carrying an annotation that
## identifies it. This is what catches a mis-assignment of labels or of
## annotations to nodes.
tagged_tree_file <- function() {
    f <- tempfile(fileext = ".tree")
    writeLines(c(
        "#NEXUS", "", "Begin trees;", "\tTranslate",
        "\t\t1 taxon_a,", "\t\t2 taxon_b,", "\t\t3 taxon_c,", "\t\t4 taxon_d",
        "\t\t;",
        paste0("\ttree TREE1 = ((3[&id=3]:1.0,1[&id=1]:1.0)[&id=101]:1.0,",
               "(4[&id=4]:1.0,2[&id=2]:1.0)[&id=102]:1.0)[&id=100]:0.0;"),
        "End;"), f)
    f
}

test_that("annotations are attached to the node they actually belong to", {
    f <- tagged_tree_file()
    tr <- read.annotated.nexus(f)
    labs <- c("taxon_a", "taxon_b", "taxon_c", "taxon_d")

    ## tip labels survive the non-TRANSLATE-order Newick
    expect_setequal(tr$tip.label, labs)
    expect_true(isTRUE(ape::all.equal.phylo(tr, ape::read.nexus(f),
                                            use.edge.length = TRUE)))

    ## each tip's annotation carries that tip's own TRANSLATE number
    for (i in seq_along(tr$tip.label)) {
        expect_equal(tr$annotations[[as.character(i)]]$id,
                     match(tr$tip.label[i], labs),
                     info = tr$tip.label[i])
    }
    ## the root, and each internal node, get their own tag
    root <- length(tr$tip.label) + 1L
    expect_equal(tr$annotations[[as.character(root)]]$id, 100)
    ca <- ape::getMRCA(tr, c("taxon_c", "taxon_a"))
    cb <- ape::getMRCA(tr, c("taxon_d", "taxon_b"))
    expect_equal(tr$annotations[[as.character(ca)]]$id, 101)
    expect_equal(tr$annotations[[as.character(cb)]]$id, 102)
})

test_that("simplify controls the shape of interval-valued annotations", {
    f <- example_tree()
    a <- read.annotated.nexus(f, simplify = TRUE)$annotations[[1]][["height_95%_HPD"]]
    b <- read.annotated.nexus(f, simplify = FALSE)$annotations[[1]][["height_95%_HPD"]]
    expect_false(is.list(a))
    expect_length(a, 2)
    expect_true(is.list(b))
    expect_equal(unlist(b), a)
})

test_that("files with many trees give a multiPhylo, no annotations and a warning", {
    f <- multi_tree_file(3)
    expect_warning(res <- read.annotated.nexus(f), "WITHOUT annotations")

    res <- suppressWarnings(read.annotated.nexus(f))
    expect_s3_class(res, "multiPhylo")
    expect_length(res, 3)
    expect_null(res[[1]]$annotations)
    expect_setequal(res[[1]]$tip.label, c("taxon_a", "taxon_b", "taxon_c"))
    expect_identical(names(res), c("STATE_1", "STATE_2", "STATE_3"))

    ## same trees as ape reads
    ref <- ape::read.nexus(f)
    expect_true(all(vapply(seq_along(ref), function(i)
        isTRUE(ape::all.equal.phylo(ref[[i]], res[[i]], use.edge.length = TRUE)),
        logical(1))))
})

test_that("a TREES block that is not closed is still read", {
    f <- multi_tree_file(3)
    x <- readLines(f)
    writeLines(x[x != "End;"], f)          # simulate a run stopped part way
    expect_warning(res <- read.annotated.nexus(f), "WITHOUT annotations")
    expect_length(res, 3)
})
