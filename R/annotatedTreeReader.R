##
## @author Marc A. Suchard
##
## A class for reading Newick formatted trees with BEAST-style annotations
##

.strip.annotations <- function(text) {
    annotations <- list()
    end <- 1

    pattern = "\\[&.*?\\]"

    repeat {
        match = regexpr(pattern=pattern,text=text)
        if (!(match[1] > 0)) {
            break
        }
        annotations[[end]] = regmatches(text, match)
        text = sub(pattern,paste("[",end,"]",sep=""), text)
        end = end + 1
    }
    return(list(annotations=annotations,tree=text))
}

.split.tree.names <- function(text) {
    text = gsub(pattern="\\[.*?\\]=",x=text,replacement="")
    text = gsub(pattern="^tree",x=text,replacement="")
    return(text)
}

.split.tree.traits <- function(text) {

    ## Pull out annotation
    text = regmatches(text,regexpr(pattern="\\[.*?\\]",text))
    ## Remove leading and trailing delimitors
    text = substring(text,3,nchar(text)-1)
    return(text)
}

.parse.value <- function(text) {
    value = text
    if (length(grep("^\\{",value))) { ## starts with {
        save = value
        value = substring(value, 2, nchar(value)-1)

        depth = 0
        r = regexpr(pattern="\\{+",value,perl=TRUE)
        match.length = attr(r, "match.length")

        if (match.length > 0) {
            depth = match.length
        }

        if (depth == 0) {
            split = ","
        } else {
            split = paste(
            "(?<=",rep("\\}",depth),")",
            ",",
            "(?=" ,rep("\\{",depth),")",
            sep="")
        }

        if (depth >= 1) {
            return(save) # TODO Still error in recursion
        }

        part = strsplit(value, split, perl=TRUE)[[1]]
        value = list()
        for (i in 1:length(part)) {
            value[[i]] = .parse.value(part[i])
        }
        ## TODO Unlist when simple array?
    } else {
        if (!is.na(suppressWarnings(as.numeric(value)))) { # is a number
            value = as.numeric(value)
        }
    }
    return(value)
}

.parse.traits <- function(text, header=FALSE) {

    if (header == TRUE) {
        text = substring(text,3,nchar(text)-1)
    }

    pattern <- "(\"[^\"]*\"+|[^,=\\s]+)\\s*(=\\s*(\\{[^=]*\\}|\"[^\"]*\"+|[^,]+))?"

    rgx <- gregexpr(pattern,text,perl=TRUE)
    n <- length(attr(rgx[[1]],"match.length"))
    traits <- list()
    start <- attr(rgx[[1]],"capture.start")
    names <- attr(rgx[[1]],"capture.names")
    length <- attr(rgx[[1]],"capture.length")
    names <- attr(rgx[[1]],"capture.names")
    for (i in 1:n) {
        s <- start[i,3]
        e <- s + length[i,3] - 1
        value <- substring(text,s,e)

        s <- start[i,1]
        e <- s + length[i,1] - 1
        key <- substring(text,s,e)

        traits[[key]] <- .parse.value(value)
    }

    return(traits)
}

.annotated.clado.build <- function(tp) {
    stop(paste("Annotated clado.build is not yet implemented.\n"))
}

## THE CODE BELOW COMES FROM 'ape'. MY GOAL IS TO DERIVE FROM THIS TO READ IN BEAST-STYLE ANNOTATIONS

.annotated.tree.build <- function(tp){

    add.internal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- current.node <<- node <<- node + 1L
        index[node] <<- j
        j <<- j + 1L
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j
        X <- unlist(strsplit(tpc[k], ":"))
        tip.label[tip] <<- X[1]
        edge.length[j] <<- as.numeric(X[2])
        
        permute[j] <<- list(.annotation.for(k))  ## permute traits

        k <<- k + 1L
        tip <<- tip + 1L
        j <<- j + 1L
    }
    go.down <- function() {
        l <- index[current.node]
        X <- unlist(strsplit(tpc[k], ":"))
        node.label[current.node - nb.tip] <<- X[1]
        edge.length[l] <<- as.numeric(X[2])
        
        permute[l] <<- list(.annotation.for(k))  ## permute traits

        k <<- k + 1L
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2L, 1L), 1, 2))
        tp <- unlist(strsplit(tp, "[\\(\\):;]"))
        obj$edge.length <- as.numeric(tp[3])
        obj$Nnode <- 1L
        obj$tip.label <- tp[2]
        if (tp[4] != "")
            obj$node.label <- tp[4]
        class(obj) <- "phylo"
        return(obj)
    }

    result = .strip.annotations(tp)
    annotations = result$annotations
    new.tp.stripped = result$tree

    annotations = lapply(annotations, .parse.traits, header=TRUE)

    tp.stripped = gsub("\\[.*?\\]","",tp)
    tpc.all <- unlist(strsplit(tp.stripped, "[\\(\\),;]"))
    keep <- nzchar(tpc.all)
    tpc <- tpc.all[keep]

    ## .strip.annotations() leaves a "[i]" placeholder wherever it took an
    ## annotation out, so the annotation belonging to the k-th node visited is
    ## found by reading the placeholder in that node's token, rather than by
    ## assuming that every node carries exactly one annotation. Trees in which
    ## only some nodes are annotated, or none at all, are handled correctly.
    ##
    ## A token is a node's own token when the delimiter that closes it is ",",
    ## ")" or ";"; the empty tokens produced by runs of "(" are not nodes. The
    ## nodes are visited in that order, so element k below belongs to visit k.
    ann.tokens <- unlist(strsplit(new.tp.stripped, "[\\(\\),;]"))
    delims <- unlist(strsplit(gsub("[^(),;]", "", tp.stripped), NULL))
    ann.index <- integer(0)
    if (length(ann.tokens) == length(tpc.all) &&
        length(delims) >= length(tpc.all)) {
        is.node <- delims[seq_along(tpc.all)] %in% c(",", ")", ";")
        node.tokens <- ann.tokens[is.node]
        idx <- suppressWarnings(as.integer(sub(".*\\[([0-9]+)\\].*", "\\1",
                                              node.tokens)))
        idx[!grepl("\\[[0-9]+\\]", node.tokens)] <- NA_integer_
        ann.index <- idx
    } else if (length(annotations)) {
        warning("could not match the annotations onto the tree: they are dropped.",
                call. = FALSE)
    }
    .annotation.for <- function(k) {
        if (k > length(ann.index) || is.na(ann.index[k])) return(NULL)
        annotations[[ann.index[k]]]
    }

    tsp <- unlist(strsplit(tp.stripped, NULL))
    skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
    nsk <- length(skeleton)
    nb.node <- sum(skeleton == ")")
    nb.tip <- sum(skeleton == ",") + 1
    nb.edge <- nb.node + nb.tip

    node.label <- character(nb.node)
    tip.label <- character(nb.tip)
    edge.length <- numeric(nb.edge)
    edge <- matrix(0L, nb.edge, 2)
    current.node <- node <- as.integer(nb.tip + 1)
    edge[nb.edge, 2] <- node
    index <- numeric(nb.edge + 1)
    index[node] <- nb.edge
    j <- k <- tip <- 1L
    
    permute = list()
    
    for (i in 2:nsk) {
        if (skeleton[i] == "(") {
            add.internal()           
        }
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") {
                add.terminal()            	
            }
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] == ",") {
                add.terminal()               
                go.down()
            }
            if (skeleton[i - 1] == ")") {
                go.down()
            }
        }        
    }
    
    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
    root.edge <- edge.length[nb.edge]
    edge.length <- edge.length[-nb.edge]
    if (!all(is.na(edge.length)))
        obj$edge.length <- edge.length
    if (is.na(node.label[1]))
        node.label[1] <- ""
    if (any(nzchar(node.label)))
        obj$node.label <- node.label
    if (!is.na(root.edge))
        obj$root.edge <- root.edge
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"

    if (length(permute) < nb.edge) permute[nb.edge] <- list(NULL)
    obj$annotations = permute
    obj
}

#' Read Newick trees with BEAST-style annotations
#'
#' Reads one or more Newick-format trees that contain BEAST-style square-
#' bracket annotations (e.g. \code{[&rate_median=0.01,length=0.5]}). The
#' annotations are parsed and stored in the \code{$annotations} list of the
#' returned \code{"phylo"} object. Largely derived from \code{ape::read.tree}.
#'
#' @param file Character. Path to the Newick file. Ignored if \code{text} is
#'   supplied. Default \code{""} (reads from standard input).
#' @param text Character vector of Newick strings, one per tree. If
#'   non-\code{NULL}, overrides \code{file}.
#' @param tree.names Character vector of tree names. Default \code{NULL}.
#' @param skip Integer. Number of lines to skip at the start of \code{file}.
#'   Default \code{0}.
#' @param comment.char Character. Lines beginning with this character are
#'   ignored. Default \code{"#"}.
#' @param keep.multi Logical. If \code{TRUE}, a single tree is still returned
#'   as a \code{"multiPhylo"} list. Default \code{FALSE}.
#' @param simplify Logical. If \code{TRUE} (default), interval-valued
#'   annotations are returned as vectors rather than as lists of length-one
#'   elements.
#' @param ... Additional arguments passed to \code{\link[base]{scan}}.
#'
#' @return A \code{"phylo"} object (or \code{"multiPhylo"} for multiple trees)
#'   with an additional \code{$annotations} list element containing the parsed
#'   BEAST annotations. As in \code{\link{read.annotated.nexus}}, the list is
#'   named by node number in the row order of \code{$edge}, with the root last.
#'
#' @seealso \code{\link{read.annotated.nexus}}, \code{\link{trann2trdat}}
#'
#' @examples
#' \dontrun{
#' tr <- read.annotated.tree("beast_mcc.tree")
#' names(tr$annotations[[1]])
#' }
#'
#' @export
read.annotated.tree <- function (file = "", text = NULL, tree.names = NULL, skip = 0,
                                  comment.char = "#", keep.multi = FALSE,
                                  simplify = TRUE, ...)
{
    unname <- function(treetext) {
        nc <- nchar(treetext)
        tstart <- 1
        while (substr(treetext, tstart, tstart) != "(" && tstart <=
               nc) tstart <- tstart + 1
        if (tstart > 1)
            return(c(substr(treetext, 1, tstart - 1), substr(treetext,
                                                             tstart, nc)))
        return(c("", treetext))
    }

    if (!is.null(text)) {
        if (!is.character(text))
            stop("argument `text' must be of mode character")
        tree <- text
    }
    else {
        tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }

    tree <- gsub("[ \n\t]", "", tree)
    tree <- gsub("\\[&R\\]", "", tree)
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    if (is.na(y[1]))
        return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree) STRING[i] <- paste(tree[x[i]:y[i]], sep = "",
                                          collapse = "")

    tmp <- unlist(lapply(STRING, unname))
    tmpnames <- tmp[c(TRUE, FALSE)]
    STRING <- tmp[c(FALSE, TRUE)]
    if (is.null(tree.names) && any(nzchar(tmpnames)))
        tree.names <- tmpnames
    colon <- grep(":", STRING)

    if (!is.null(tree.names)) {
        traits.text = lapply(tree.names, .split.tree.traits)
        tree.names = lapply(tree.names, .split.tree.names)
        tree.traits = lapply(traits.text, .parse.traits)
    }

    if (!length(colon)) {
        stop(paste("Annotated clado.build is not yet implemented.\n"))
        obj <- lapply(STRING, .annotated.clado.build)
    }
    else if (length(colon) == Ntree) {
        obj <- lapply(STRING, .annotated.tree.build)
    }
    else {
        obj <- vector("list", Ntree)
        obj[colon] <- lapply(STRING[colon], .annotated.tree.build)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        obj[nocolon] <- lapply(STRING[nocolon], clado.build)
    }
    for (i in 1:Ntree) {
        ROOT <- length(obj[[i]]$tip.label) + 1
        if (sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] >
            1)
            stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading Newick file aborted at tree no.",
                       i))
    }
    ## name the annotations by the node they describe, as read.annotated.nexus does
    for (i in 1:Ntree) {
        if (!is.null(obj[[i]]$annotations)) {
            obj[[i]] <- .name.annotations(obj[[i]])
            if (simplify) {
                obj[[i]]$annotations <- .simplify.annotations(obj[[i]]$annotations)
            }
        }
    }
    if (Ntree == 1 && !keep.multi)
        obj <- obj[[1]]
    else {
        if (!is.null(tree.names)) {
            names(obj) <- tree.names
        }
        class(obj) <- "multiPhylo"
    }
    obj
}

## Order in which the tips occur in a Newick string.
##
## .annotated.tree.build() numbers the tips in the order it meets them while
## walking the string, so this is the order of the tip indices used in $edge.
## A tip label is a token that follows "(" or "," and is followed by ":"; an
## internal node label follows ")" and is therefore not matched.
.tip.order <- function(tp) {
    s <- gsub("\\[.*?\\]", "", tp)
    s <- gsub("[[:space:]]", "", s)
    hits <- regmatches(s, gregexpr("[(,][^(),:]+:", s))[[1]]
    sub(":$", "", sub("^[(,]", "", hits))
}


## Name the elements of tree$annotations with the node they describe, following
## the row order of tree$edge. NELSI appends the root's annotations last.
.name.annotations <- function(tree) {
    n.edge <- nrow(tree$edge)
    n.ann <- length(tree$annotations)
    root <- length(tree$tip.label) + 1L
    if (n.ann == n.edge + 1L) {
        nodes <- c(tree$edge[, 2], root)
    } else if (n.ann == n.edge) {
        nodes <- tree$edge[, 2]
    } else {
        stop(sprintf("cannot match %d annotations to %d edges.\n", n.ann, n.edge))
    }
    names(tree$annotations) <- as.character(nodes)
    tree
}


## Turn the two-element lists that .parse.traits returns for interval-valued
## annotations (height_95%_HPD, length_range, ...) into plain vectors.
.simplify.annotations <- function(annotations) {
    lapply(annotations, function(a) {
        if (!is.list(a)) return(a)
        lapply(a, function(v) if (is.list(v) && all(lengths(v) == 1L)) unlist(v) else v)
    })
}


#' Read a NEXUS file with BEAST-style annotations
#'
#' Reads a NEXUS-format file produced by BEAST or similar software, parsing
#' branch-level annotations stored in BEAST's square-bracket notation and
#' handling TRANSLATE blocks.
#'
#' The annotations are returned in \code{$annotations}, a \emph{named} list
#' whose names are the node numbers in the row order of \code{$edge}:
#'
#' \preformatted{names(tree$annotations)[i] == as.character(tree$edge[i, 2])}
#'
#' Element \code{i} therefore holds the annotations of the node at the child end
#' of edge \code{i}, describing both that node and the branch
#' \code{tree$edge.length[i]} subtending it. One extra, final element holds the
#' annotations of the root and is named with the root node number
#' (\code{Ntip(tree) + 1}). Elements can be retrieved either positionally
#' (\code{tree$annotations[[i]]}, matching edge \code{i}) or by node number
#' (\code{tree$annotations[["36"]]}).
#'
#' Annotations are only parsed for files that contain a single tree, which is
#' the usual case for summary trees such as those written by TreeAnnotator. If
#' the file holds several trees a \code{"multiPhylo"} object is returned
#' \emph{without} annotations and a warning is issued.
#'
#' @section Tip labels:
#' Earlier versions assigned translated tip labels positionally
#' (\code{tree$tip.label <- TRANS[, 2]}), which assumes the tips appear in the
#' Newick string in the same order as in the TRANSLATE block. BEAST writes them
#' in topological order instead, so those versions returned a tree with the
#' correct topology and branch lengths but permuted, i.e. wrong, tip labels.
#' The order in which the tips actually occur in the Newick string is now
#' re-derived and the labels are mapped through the TRANSLATE table.
#'
#' @param file Character. Path to the NEXUS file.
#' @param tree.names Character vector of tree names to assign. Default
#'   \code{NULL} (names are read from the file).
#' @param simplify Logical. If \code{TRUE} (default), interval-valued
#'   annotations such as \code{height_95\%_HPD} are returned as length-two
#'   vectors rather than as lists of two length-one elements.
#' @param check Logical. If \code{TRUE}, the tree is compared with the one
#'   returned by \code{ape::read.nexus} and a warning is issued if they differ.
#'   Default \code{FALSE}.
#'
#' @return A \code{"phylo"} object with an additional named \code{$annotations}
#'   list, as described above. For files containing more than one tree, a
#'   \code{"multiPhylo"} object with no annotations.
#'
#' @seealso \code{\link{read.annotated.tree}}, \code{\link{trann2trdat}}
#'
#' @examples
#' \dontrun{
#' tr <- read.annotated.nexus("beast_mcc.tree")
#' tr$annotations[[1]]                       # annotations of node tr$edge[1, 2]
#' tr$annotations[[as.character(Ntip(tr) + 1)]]   # annotations of the root
#' dat <- trann2trdat(tr)
#' }
#'
#' @export
read.annotated.nexus <- function (file, tree.names = NULL, simplify = TRUE,
                                  check = FALSE) {
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)

    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    if (!length(i1)) stop("no TREES block found in ", file, "\n")
    i1 <- i1[1]
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- length(i2) == 1 && i2 > i1
    if (translation) {
        end <- semico[semico > i2][1]
        x <- unlist(strsplit(X[(i2 + 1):end], "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    start <- if (translation) semico[semico > i2][1] + 1 else semico[semico > i1][1]
    if (is.na(start)) stop("could not find the start of the TREES block in ", file, "\n")
    ## a file whose TREES block is not closed (e.g. a run stopped part way
    ## through writing) is read up to its last line rather than failing
    end <- endblock[endblock > i1][1] - 1
    if (is.na(end)) end <- length(X)
    tree <- X[start:end]
    rm(X)

    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico)
    if (Ntree == 1 && length(tree) > 1) {
        STRING <- paste(tree, collapse = "")
    } else {
        if (any(diff(semico) != 1)) {
            STRING <- character(Ntree)
            s <- c(1, semico[-Ntree] + 1)
            j <- mapply(":", s, semico)
            if (is.list(j)) {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], collapse = "")
            } else {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[, i]], collapse = "")
            }
        } else STRING <- tree
    }
    rm(tree)

    STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
    Ntree <- length(STRING)
    if (Ntree == 0) stop("no tree found in the TREES block.\n")

    ## Annotations are only parsed for single-tree files. For tree sets (e.g. a
    ## posterior sample) fall back on ape::read.nexus, which reads the topology
    ## and branch lengths but discards the [&...] metadata.
    nms.trees <- sub("[[:blank:]]*=.*", "", STRING)
    nms.trees <- sub("^[[:blank:]]*tree[[:blank:]]*", "", nms.trees, ignore.case = TRUE)

    STRING <- gsub("\\[&R\\]", "", STRING)
    STRING <- sub("^.*?= *", "", STRING)
    STRING <- gsub("\\s", "", STRING)

    if (Ntree > 1) {
        warning(sprintf(paste0("the file contains %d trees: returning a ",
                               "'multiPhylo' object WITHOUT annotations. ",
                               "Annotations are only parsed for files with a ",
                               "single tree, such as summary trees written by ",
                               "TreeAnnotator."), Ntree), call. = FALSE)
        ## strip the annotations and let ape's parser build the trees
        tp <- gsub("\\[.*?\\]", "", STRING)
        keep <- grepl(";$", tp)
        if (any(!keep)) {
            warning(sprintf("dropping %d incomplete tree(s) at the end of the file.",
                            sum(!keep)), call. = FALSE)
        }
        if (!any(keep)) stop("no complete tree found in ", file, "\n")
        trees <- ape::read.tree(text = tp[keep])
        if (inherits(trees, "phylo")) {
            trees <- structure(list(trees), class = "multiPhylo")
        }
        if (translation) {
            for (i in seq_along(trees)) {
                ind <- match(trees[[i]]$tip.label, TRANS[, 1])
                if (anyNA(ind)) {
                    stop("tip(s) missing from the TRANSLATE block: ",
                         paste(trees[[i]]$tip.label[is.na(ind)], collapse = ", "), "\n")
                }
                trees[[i]]$tip.label <- TRANS[ind, 2]
            }
        }
        ## use ape's shared tip-label convention when the trees share a tip set
        cmp <- try(ape::.compressTipLabel(trees), silent = TRUE)
        if (!inherits(cmp, "try-error")) trees <- cmp
        if (!is.null(tree.names)) {
            names(trees) <- tree.names
        } else if (!all(nms.trees[keep] == "")) {
            names(trees) <- nms.trees[keep]
        }
        return(trees)
    }

    if (!length(grep(":", STRING))) {
        stop(".annotated.clado.build is not yet implemented.\n")
    }

    trees <- .annotated.tree.build(STRING)

    if (!translation) n <- length(trees$tip.label)
    ROOT <- n + 1
    if (sum(trees$edge[, 1] == ROOT) == 1 && dim(trees$edge)[1] > 1) {
        stop("The tree has apparently singleton node(s): cannot read tree file.\n")
    }

    ## Map the tips onto their labels using the order in which they occur in the
    ## Newick string, which is the order .annotated.tree.build numbered them in.
    if (translation) {
        tips <- .tip.order(STRING)
        if (length(tips) != length(trees$tip.label)) {
            stop(sprintf(paste0("found %d tips in the Newick string but the ",
                                "parsed tree has %d: cannot assign tip labels.\n"),
                         length(tips), length(trees$tip.label)))
        }
        ind <- match(tips, TRANS[, 1])
        if (anyNA(ind)) {
            stop("tip(s) missing from the TRANSLATE block: ",
                 paste(tips[is.na(ind)], collapse = ", "), "\n")
        }
        trees$tip.label <- TRANS[ind, 2]
    }

    if (!is.null(trees$annotations)) {
        trees <- .name.annotations(trees)
        if (simplify) trees$annotations <- .simplify.annotations(trees$annotations)
    }

    if (check) {
        ref <- try(ape::read.nexus(file), silent = TRUE)
        if (!inherits(ref, "try-error") &&
            !isTRUE(ape::all.equal.phylo(trees, ref, use.edge.length = FALSE))) {
            warning("the tree does not match the one returned by ape::read.nexus(); ",
                    "check the file.", call. = FALSE)
        }
    }

    trees
} # end read.annotated.nexus
