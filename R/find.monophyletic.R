#' Find monophyletic clades matching a tip-label pattern
#'
#' Identifies all maximal monophyletic clades whose tip labels all match a
#' given pattern. The algorithm traverses the tree from the root downward,
#' pruning subtrees that have no matching tips, and collecting subtrees where
#' every tip matches.
#'
#' @param tr A rooted phylogenetic tree of class \code{"phylo"}.
#' @param tag Character string (regular expression) to match against tip
#'   labels.
#' @param include.singletons Logical. If \code{TRUE}, tips that match
#'   \code{tag} but do not belong to any monophyletic multi-tip clade are also
#'   included in the output. Default \code{FALSE}.
#'
#' @return A list of character vectors. Each element is a vector of tip labels
#'   belonging to one monophyletic clade. If \code{include.singletons = TRUE},
#'   singleton matching tips are appended as individual character scalars.
#'
#' @seealso \code{\link{get.mrca}}, \code{\link{is.polytomy}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rtree(12)
#' tr$tip.label[c(1, 2, 5)] <- paste0("TAG_", tr$tip.label[c(1, 2, 5)])
#' find.monophyletic(tr, "TAG_")
#' find.monophyletic(tr, "TAG_", include.singletons = TRUE)
#'
#' @export
find.monophyletic <- function(tr, tag, include.singletons = FALSE) {
    if (!ape::is.rooted(tr)) stop('tree is not rooted')
    Ntips <- length(tr$tip.label)
    tips <- seq_len(Ntips)
    intnodes <- unique(tr$edge[, 1])

    # Replaces geiger's internal .Call("get_descendants") with castor
    get_descendants <- function(ape_node) {
        sub <- castor::get_subtree_at_node(tr, ape_node - Ntips)
        c(sub$new2old_tip, sub$new2old_node + Ntips)
    }

    clades <- list()
    i <- 1L
    while (i <= length(intnodes)) {
        descending_nodes <- get_descendants(intnodes[i])
        descending_tips  <- tr$tip.label[descending_nodes[descending_nodes %in% tips]]
        num_matches      <- grep(tag, descending_tips)

        if (length(num_matches) == 0L) {
            intnodes <- intnodes[!intnodes %in% descending_nodes]
        } else if (length(num_matches) == length(descending_tips)) {
            clades[[length(clades) + 1L]] <- descending_tips
            intnodes <- intnodes[!intnodes %in% descending_nodes]
        } else if (is.polytomy(tr, intnodes[i])) {
            i_matches       <- which(tr$tip.label %in% descending_tips[num_matches])
            ancestral_nodes <- sapply(i_matches, function(x) get.mrca(tr, x))
            tips_in_multi_clade <- i_matches[ancestral_nodes == intnodes[i]]
            if (length(tips_in_multi_clade) > 1L) {
                clades[[length(clades) + 1L]] <- tr$tip.label[tips_in_multi_clade]
            }
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
