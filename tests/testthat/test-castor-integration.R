library(testthat)
library(ape)
library(castor)
library(phangorn)

# Source all NELSI R files so tests run without installing the package
nelsi_r <- list.files(
    system.file("../../../R", package = "NELSI", mustWork = FALSE),
    full.names = TRUE, pattern = "\\.R$"
)
if (length(nelsi_r) == 0) {
    # Running from repo root during development
    nelsi_r <- list.files(
        file.path(dirname(dirname(normalizePath(testthat::test_path()))), "R"),
        full.names = TRUE, pattern = "\\.R$"
    )
}
invisible(lapply(nelsi_r, source))

# ---------------------------------------------------------------------------
# Old (pre-castor) reference implementations defined inline
# ---------------------------------------------------------------------------

.allnode_times_old <- function(phylo, tipsonly = FALSE, reverse = TRUE, keeproot = FALSE) {
    di.phylo   <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times  <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, ]
    if (reverse) {
        node.times <- abs(node.times - max(node.times))
        if (keeproot) node.times <- node.times + phylo$root.edge
    }
    if (tipsonly) node.times <- node.times[names(node.times) %in% 1:length(phylo$tip.label)]
    return(node.times)
}

.intnode_times_old <- function(phylo) {
    di.phylo   <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times  <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo,
                                           (length(phylo$tip.label) + 1):nrow(di.phylo)]
    return(node.times)
}

.get_rtt_dist_old <- function(phylo) {
    di.phylo   <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times  <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, ]
    tip.times   <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo,
                                           1:length(phylo$tip.label)]
    tip.times   <- abs(tip.times - max(node.times))
    return(tip.times)
}

.node_to_tip_dist_old <- function(phylo) {
    di.phylo   <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times  <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo,
                                           1:length(phylo$tip.label)]
    node.to.tip.distance <- abs(node.times - max(.intnode_times_old(phylo)))
    return(node.to.tip.distance)
}

.get_mrca_old <- function(tr, tips) {
    root <- tr$edge[, 1][!tr$edge[, 1] %in% tr$edge[, 2]][1]
    paths_to_root <- list()
    for (i in seq_along(tips)) {
        nodes_in_path <- tips[i]
        ancestor <- 0L
        while (ancestor != root) {
            ancestor <- tr$edge[nodes_in_path[length(nodes_in_path)] == tr$edge[, 2], 1]
            nodes_in_path <- c(nodes_in_path, ancestor)
        }
        paths_to_root[[i]] <- nodes_in_path
    }
    rintersect <- function(x) {
        if (length(x) == 2L) return(intersect(x[[1]], x[[2]]))
        return(intersect(x[[1]], rintersect(x[-1])))
    }
    if (length(tips) > 1L) return(rintersect(paths_to_root)[1])
    return(paths_to_root[[1]][2])
}

.get_ancestor_old <- function(tr, target_node) {
    root_node <- unique(tr$edge[!(tr$edge[, 1] %in% tr$edge[, 2]), 1])
    if (length(target_node) > 1L) target_node <- getMRCA(tr, target_node)
    ancestors_nodes    <- target_node
    ancestors_branches <- vector()
    while (!(root_node %in% ancestors_nodes)) {
        matches <- tr$edge[, 2] == ancestors_nodes[length(ancestors_nodes)]
        ancestors_nodes    <- c(ancestors_nodes, tr$edge[matches, 1])
        ancestors_branches <- c(ancestors_branches, which(matches))
    }
    list(ancestor.nodes = ancestors_nodes, ancestor.branches = ancestors_branches)
}

.get_descending_old <- function(tr, target_node) {
    all_tips <- tr$edge[!(tr$edge[, 2] %in% tr$edge[, 1]), 2]
    descendant_nodes    <- target_node
    descendant_branches <- vector()
    nodes_temp <- target_node
    while (!all(nodes_temp %in% all_tips)) {
        matches    <- tr$edge[, 1] %in% nodes_temp
        nodes_temp <- tr$edge[matches, 2]
        descendant_nodes    <- c(descendant_nodes, nodes_temp)
        descendant_branches <- c(descendant_branches, which(matches))
    }
    list(descending.nodes = descendant_nodes, descending.branches = descendant_branches)
}

.tips_in_clade <- function(tr, node) {
    # Replicates phangorn::tips() which is no longer exported in current phangorn
    Ntips <- length(tr$tip.label)
    desc <- node; temp <- node
    repeat {
        ch <- tr$edge[tr$edge[, 1] %in% temp, 2]
        if (length(ch) == 0L) break
        desc <- c(desc, ch); temp <- ch  # accumulate ALL descendants
    }
    tr$tip.label[desc[desc <= Ntips]]
}

.sum_descending_old <- function(tr, target_node) {
    # Replicates old sum.descending.branches without phangorn functions
    # (phangorn::tips and phangorn::nodepath removed from exports in current versions)
    Ntips   <- length(tr$tip.label)
    parent  <- integer(Ntips + tr$Nnode)
    parent[tr$edge[, 2]] <- tr$edge[, 1]

    tip_labels       <- .tips_in_clade(tr, target_node)
    tips_descendants <- match(tip_labels, tr$tip.label)

    path_to_node <- function(tip) {
        path <- tip
        cur  <- tip
        while (cur != target_node) { cur <- parent[cur]; path <- c(path, cur) }
        path
    }

    descending_nodes <- unique(unlist(lapply(tips_descendants, path_to_node)))
    descending_nodes <- descending_nodes[descending_nodes != target_node]
    branch_lengths   <- tr$edge.length[tr$edge[, 2] %in% descending_nodes]
    sum(branch_lengths)
}

.path.node_old <- function(phylo, tipsonly = TRUE) {
    di.tr   <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    if (tipsonly) {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == root.tr, 1:length(phylo$tip.label)]
        nodesinpath   <- sapply(seq_len(length(phylo$tip.label)),
                                function(x) length(phangorn::Ancestors(phylo, x)))
    } else {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == root.tr, ]
        nodesinpath   <- sapply(seq_len(length(phylo$tip.label) + phylo$Nnode),
                                function(x) length(phangorn::Ancestors(phylo, x)))
    }
    list(roottotippath = roottotippath, nodesinpath = nodesinpath)
}

# ---------------------------------------------------------------------------
# Test trees (fixed seeds for reproducibility)
# ---------------------------------------------------------------------------

set.seed(42)
tr_ultra    <- rcoal(10)       # ultrametric 10-tip
tr_nonultra <- rtree(10)       # non-ultrametric 10-tip
tr_large    <- rcoal(50)       # larger ultrametric tree

set.seed(42)
tr_large_nonultra <- rtree(50) # larger non-ultrametric tree

# Tree with tags for find.monophyletic
set.seed(99)
tr_tagged <- rtree(14)
tr_tagged$tip.label[c(1, 2, 4, 9, 11)] <- paste0("TAG_", tr_tagged$tip.label[c(1, 2, 4, 9, 11)])

# ---------------------------------------------------------------------------
# Step 2 — all.node.times
# ---------------------------------------------------------------------------

test_that("all.node.times matches old output on ultrametric tree", {
    old <- .allnode_times_old(tr_ultra)
    new <- all.node.times(tr_ultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("all.node.times matches old output on non-ultrametric tree", {
    old <- .allnode_times_old(tr_nonultra)
    new <- all.node.times(tr_nonultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("all.node.times tipsonly=TRUE returns only tips", {
    new <- all.node.times(tr_ultra, tipsonly = TRUE)
    expect_length(new, length(tr_ultra$tip.label))
    old <- .allnode_times_old(tr_ultra, tipsonly = TRUE)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("all.node.times reverse=FALSE matches old output", {
    old <- .allnode_times_old(tr_nonultra, reverse = FALSE)
    new <- all.node.times(tr_nonultra, reverse = FALSE)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("all.node.times scales to large tree", {
    old <- .allnode_times_old(tr_large)
    new <- all.node.times(tr_large)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Step 4 — int.node.times
# ---------------------------------------------------------------------------

test_that("int.node.times matches old output on ultrametric tree", {
    old <- .intnode_times_old(tr_ultra)
    new <- int.node.times(tr_ultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("int.node.times matches old output on non-ultrametric tree", {
    old <- .intnode_times_old(tr_nonultra)
    new <- int.node.times(tr_nonultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("int.node.times returns Nnodes values", {
    new <- int.node.times(tr_ultra)
    expect_length(new, tr_ultra$Nnode)
})

# ---------------------------------------------------------------------------
# Step 3 — get.rtt.dist and node.to.tip.dist
# ---------------------------------------------------------------------------

test_that("get.rtt.dist matches old output on ultrametric tree", {
    old <- .get_rtt_dist_old(tr_ultra)
    new <- get.rtt.dist(tr_ultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("get.rtt.dist matches old output on non-ultrametric tree", {
    old <- .get_rtt_dist_old(tr_nonultra)
    new <- get.rtt.dist(tr_nonultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("node.to.tip.dist matches old output on ultrametric tree", {
    old <- .node_to_tip_dist_old(tr_ultra)
    new <- node.to.tip.dist(tr_ultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("node.to.tip.dist matches old output on non-ultrametric tree", {
    old <- .node_to_tip_dist_old(tr_nonultra)
    new <- node.to.tip.dist(tr_nonultra)
    expect_equal(as.numeric(old), as.numeric(new), tolerance = 1e-10)
})

test_that("get.rtt.dist and node.to.tip.dist agree with each other", {
    expect_equal(as.numeric(get.rtt.dist(tr_nonultra)),
                 as.numeric(node.to.tip.dist(tr_nonultra)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Step 5 — get.mrca
# ---------------------------------------------------------------------------

test_that("get.mrca matches old for two tips (ultrametric)", {
    for (pair in list(c(1,2), c(3,7), c(1,10), c(5,6))) {
        old <- .get_mrca_old(tr_ultra, pair)
        new <- get.mrca(tr_ultra, pair)
        expect_equal(new, old, label = paste("pair:", paste(pair, collapse=",")))
    }
})

test_that("get.mrca matches old for two tips (non-ultrametric)", {
    for (pair in list(c(1,2), c(4,8), c(2,9))) {
        old <- .get_mrca_old(tr_nonultra, pair)
        new <- get.mrca(tr_nonultra, pair)
        expect_equal(new, old, label = paste("pair:", paste(pair, collapse=",")))
    }
})

test_that("get.mrca single tip returns parent node", {
    for (tip in 1:5) {
        old <- .get_mrca_old(tr_ultra, tip)
        new <- get.mrca(tr_ultra, tip)
        expect_equal(new, old, label = paste("tip:", tip))
    }
})

test_that("get.mrca matches old for multi-tip set", {
    tips <- c(1, 3, 5, 8)
    old  <- .get_mrca_old(tr_ultra, tips)
    new  <- get.mrca(tr_ultra, tips)
    expect_equal(new, old)
})

# ---------------------------------------------------------------------------
# Step 6 — get.ancestor.nodes.branches
# ---------------------------------------------------------------------------

test_that("get.ancestor.nodes.branches: node set matches old (order-insensitive)", {
    Ntips  <- length(tr_ultra$tip.label)
    target <- Ntips + 3L   # an internal node
    old <- .get_ancestor_old(tr_ultra, target)
    new <- get.ancestor.nodes.branches(tr_ultra, target)
    expect_setequal(new$ancestor.nodes, old$ancestor.nodes)
    expect_setequal(new$ancestor.branches, old$ancestor.branches)
})

test_that("get.ancestor.nodes.branches: root is always last node", {
    Ntips    <- length(tr_ultra$tip.label)
    root_node <- tr_ultra$edge[!(tr_ultra$edge[,1] %in% tr_ultra$edge[,2]), 1][1]
    for (internal in seq(Ntips + 1L, Ntips + tr_ultra$Nnode)) {
        res <- get.ancestor.nodes.branches(tr_ultra, internal)
        expect_equal(res$ancestor.nodes[length(res$ancestor.nodes)], root_node)
    }
})

test_that("get.ancestor.nodes.branches: target_node is first in ancestor.nodes", {
    Ntips  <- length(tr_nonultra$tip.label)
    target <- Ntips + 4L
    res    <- get.ancestor.nodes.branches(tr_nonultra, target)
    expect_equal(res$ancestor.nodes[1], target)
})

test_that("get.ancestor.nodes.branches: branch count = node count - 1", {
    Ntips  <- length(tr_ultra$tip.label)
    target <- Ntips + 2L
    res    <- get.ancestor.nodes.branches(tr_ultra, target)
    expect_equal(length(res$ancestor.branches), length(res$ancestor.nodes) - 1L)
})

test_that("get.ancestor.nodes.branches: branches index correct edges", {
    Ntips  <- length(tr_ultra$tip.label)
    target <- Ntips + 5L
    res    <- get.ancestor.nodes.branches(tr_ultra, target)
    # Each branch index should point to an edge whose child = corresponding ancestor node
    # The i-th branch leads from ancestor[i+1] to ancestor[i]
    for (k in seq_along(res$ancestor.branches)) {
        edge_child <- tr_ultra$edge[res$ancestor.branches[k], 2]
        expect_equal(edge_child, res$ancestor.nodes[k])
    }
})

# ---------------------------------------------------------------------------
# Step 7 — get.descending.nodes.branches and sum.descending.branches
# ---------------------------------------------------------------------------

test_that("get.descending.nodes.branches: node set matches old", {
    Ntips  <- length(tr_ultra$tip.label)
    target <- Ntips + 2L
    old <- .get_descending_old(tr_ultra, target)
    new <- get.descending.nodes.branches(tr_ultra, target)
    expect_setequal(new$descending.nodes, old$descending.nodes)
    expect_setequal(new$descending.branches, old$descending.branches)
})

test_that("get.descending.nodes.branches: all tips within set are descendants", {
    Ntips  <- length(tr_nonultra$tip.label)
    for (internal in seq(Ntips + 1L, Ntips + tr_nonultra$Nnode)) {
        res   <- get.descending.nodes.branches(tr_nonultra, internal)
        n_sub <- castor::get_subtree_at_node(tr_nonultra, internal - Ntips)
        expect_setequal(
            sort(n_sub$new2old_tip),
            sort(res$descending.nodes[res$descending.nodes <= Ntips])
        )
    }
})

test_that("get.descending.nodes.branches: tip node returns itself, no branches", {
    res <- get.descending.nodes.branches(tr_ultra, 1L)
    expect_equal(res$descending.nodes, 1L)
    expect_length(res$descending.branches, 0L)
})

test_that("sum.descending.branches matches old on ultrametric tree", {
    Ntips <- length(tr_ultra$tip.label)
    for (internal in seq(Ntips + 1L, Ntips + tr_ultra$Nnode)) {
        old <- .sum_descending_old(tr_ultra, internal)
        new <- sum.descending.branches(tr_ultra, internal)
        expect_equal(new, old, tolerance = 1e-10,
                     label = paste("internal node:", internal))
    }
})

test_that("sum.descending.branches matches old on non-ultrametric tree", {
    Ntips <- length(tr_nonultra$tip.label)
    for (internal in seq(Ntips + 1L, Ntips + tr_nonultra$Nnode)) {
        old <- .sum_descending_old(tr_nonultra, internal)
        new <- sum.descending.branches(tr_nonultra, internal)
        expect_equal(new, old, tolerance = 1e-10,
                     label = paste("internal node:", internal))
    }
})

test_that("sum.descending.branches: tip node returns 0", {
    expect_equal(sum.descending.branches(tr_ultra, 1L), 0)
})

# ---------------------------------------------------------------------------
# Step 8 — get.ltt.summary (unchanged; regression test only)
# ---------------------------------------------------------------------------

test_that("get.ltt.summary errors on ultrametric tree", {
    expect_error(get.ltt.summary(tr_ultra), "ultrametric")
})

test_that("get.ltt.summary returns named vector of length 4 on non-ultrametric tree", {
    res <- get.ltt.summary(tr_nonultra)
    expect_length(res, 4L)
    expect_named(res, c("time_max_lineages", "slope1", "slope2", "ratio_slopes"))
})

test_that("get.ltt.summary: ratio_slopes == slope1/slope2", {
    res <- get.ltt.summary(tr_nonultra)
    expect_equal(unname(res["ratio_slopes"]),
                 unname(res["slope1"] / res["slope2"]), tolerance = 1e-10)
})

test_that("get.ltt.summary: relative_time_max_l is in (0, 1]", {
    res <- get.ltt.summary(tr_large_nonultra)
    expect_true(res["time_max_lineages"] > 0 && res["time_max_lineages"] <= 1)
})

# ---------------------------------------------------------------------------
# Step 9 — find.monophyletic
# ---------------------------------------------------------------------------

# Helper: normalise a list of tip-label vectors for order-insensitive comparison
.norm_clades <- function(clades) {
    sorted <- lapply(clades, sort)
    sorted[order(sapply(sorted, function(x) x[1]))]
}

test_that("find.monophyletic returns same clades as geiger-based original", {
    skip_if_not_installed("geiger")
    old <- suppressWarnings(suppressMessages({
        geiger_find_mono <- find.monophyletic  # will be overwritten; use inline old version
        # Inline old implementation
        local({
            find.monophyletic_old_impl <- function(tr, tag, include.singletons = FALSE) {
                require(geiger, quietly = TRUE)
                get.descendants.of.node <- function(phy, node, tips = FALSE) {
                    n   <- Ntip(phy)
                    all <- ifelse(tips, FALSE, TRUE)
                    out <- .Call("get_descendants", tree = list(
                        NODE = as.integer(node), ROOT = as.integer(n + 1),
                        ALL  = as.integer(all),
                        ENDOFCLADE = as.integer(dim(phy$edge)[1]),
                        ANC = as.integer(phy$edge[, 1]),
                        DES = as.integer(phy$edge[, 2])), PACKAGE = "geiger")
                    res <- out$TIPS
                    if (!length(res)) res <- NULL
                    return(c(node, res))
                }
                if (!is.rooted(tr)) stop('tree is not rooted')
                tips_idx <- 1:length(tr$tip.label)
                intnodes <- unique(tr$edge[, 1])
                clades   <- list()
                i <- 1L
                while (i <= length(intnodes)) {
                    descending_nodes <- get.descendants.of.node(tr, intnodes[i])
                    descending_tips  <- tr$tip.label[descending_nodes[descending_nodes %in% tips_idx]]
                    num_matches <- grep(tag, descending_tips)
                    if (length(num_matches) == 0L) {
                        intnodes <- intnodes[-which(intnodes %in% descending_nodes)]
                    } else if (length(num_matches) == length(descending_tips)) {
                        clades[[length(clades) + 1L]] <- descending_tips
                        intnodes <- intnodes[-which(intnodes %in% descending_nodes)]
                    } else if ({d <- dim(tr$edge[tr$edge[,1]==intnodes[i],]); d[1] > 2}) {
                        i_matches       <- which(tr$tip.label %in% descending_tips[num_matches])
                        ancestral_nodes <- sapply(i_matches, function(x) .get_mrca_old(tr, x))
                        tips_in_multi_clade <- i_matches[ancestral_nodes == intnodes[i]]
                        if (length(tips_in_multi_clade) > 1L)
                            clades[[length(clades) + 1L]] <- tr$tip.label[tips_in_multi_clade]
                        intnodes <- intnodes[intnodes != intnodes[i]]
                    } else {
                        i <- i + 1L
                    }
                }
                if (include.singletons) {
                    tips_in_clades  <- unlist(clades)
                    tips_out_clades <- tr$tip.label[!tr$tip.label %in% tips_in_clades]
                    tips_out_clades <- grep(tag, tips_out_clades, value = TRUE)
                    return(c(clades, tips_out_clades))
                }
                return(clades)
            }
            find.monophyletic_old_impl(tr_tagged, "TAG_")
        })
    }))
    new <- find.monophyletic(tr_tagged, "TAG_")
    expect_equal(.norm_clades(old), .norm_clades(new))
})

test_that("find.monophyletic: all returned tip labels match the tag", {
    res <- find.monophyletic(tr_tagged, "TAG_")
    all_returned <- unlist(res)
    expect_true(all(grepl("TAG_", all_returned)))
})

test_that("find.monophyletic: singletons included when requested", {
    with_sing    <- find.monophyletic(tr_tagged, "TAG_", include.singletons = TRUE)
    without_sing <- find.monophyletic(tr_tagged, "TAG_", include.singletons = FALSE)
    # With singletons has >= as many tip labels
    expect_gte(length(unlist(with_sing)), length(unlist(without_sing)))
    # All returned labels match tag
    expect_true(all(grepl("TAG_", unlist(with_sing))))
})

test_that("find.monophyletic: errors on unrooted tree", {
    tr_unrooted <- unroot(tr_ultra)
    expect_error(find.monophyletic(tr_unrooted, "t"), "rooted")
})

test_that("find.monophyletic: returns empty list when no tip matches tag", {
    res <- find.monophyletic(tr_ultra, "NOTHERE_XYZ")
    expect_length(res, 0L)
})

# ---------------------------------------------------------------------------
# Step 10 — path.node
# ---------------------------------------------------------------------------

test_that("path.node tipsonly=TRUE: root-to-tip distances match old", {
    old <- .path.node_old(tr_ultra, tipsonly = TRUE)
    pdf(NULL)
    new <- path.node(tr_ultra, tipsonly = TRUE)
    dev.off()
    expect_equal(as.numeric(old$roottotippath), as.numeric(new$roottotippath),
                 tolerance = 1e-10)
})

test_that("path.node tipsonly=TRUE: ancestor counts match old", {
    old <- .path.node_old(tr_ultra, tipsonly = TRUE)
    pdf(NULL)
    new <- path.node(tr_ultra, tipsonly = TRUE)
    dev.off()
    expect_equal(as.numeric(old$nodesinpath), as.numeric(new$nodesinpath))
})

test_that("path.node tipsonly=FALSE: distances match old (all nodes)", {
    old <- .path.node_old(tr_ultra, tipsonly = FALSE)
    pdf(NULL)
    new <- path.node(tr_ultra, tipsonly = FALSE)
    dev.off()
    expect_equal(as.numeric(old$roottotippath), as.numeric(new$roottotippath),
                 tolerance = 1e-10)
})

test_that("path.node tipsonly=FALSE: ancestor counts match old (all nodes)", {
    old <- .path.node_old(tr_nonultra, tipsonly = FALSE)
    pdf(NULL)
    new <- path.node(tr_nonultra, tipsonly = FALSE)
    dev.off()
    expect_equal(as.numeric(old$nodesinpath), as.numeric(new$nodesinpath))
})

test_that("path.node: root has 0 ancestors", {
    Ntips     <- length(tr_ultra$tip.label)
    root_node <- tr_ultra$edge[!(tr_ultra$edge[,1] %in% tr_ultra$edge[,2]), 1][1]
    pdf(NULL)
    res <- path.node(tr_ultra, tipsonly = FALSE)
    dev.off()
    expect_equal(res$nodesinpath[root_node], 0L)
})

# ---------------------------------------------------------------------------
# Integration smoke test — simulate.rate still works end-to-end
# ---------------------------------------------------------------------------

test_that("simulate.rate end-to-end run completes without error", {
    # simulate.rate(tree, FUN, ...) — second arg is the simulation function itself
    set.seed(1)
    tr <- rcoal(10)
    expect_no_error({
        rate_sim <- simulate.rate(tr, simulate.autocor.kishino,
                                  list(initial.rate = 0.01, v = 0.3))
    })
    expect_s3_class(rate_sim, "ratesim")
    expect_named(rate_sim, c("phylogram", "tree.data.matrix"))
})
