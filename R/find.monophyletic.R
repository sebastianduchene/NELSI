find.monophyletic <- function(tr, tag, include.singletons = F){
    require(geiger)
    get.descendants.of.node <-
        function(phy, node, tips=FALSE){
            n=Ntip(phy)
            all=ifelse(tips, FALSE, TRUE)
            out <- .Call("get_descendants", tree=list(
                                                NODE = as.integer(node),
                                                ROOT = as.integer(n+1),
                                                ALL = as.integer(all),
                                                ENDOFCLADE = as.integer(dim(phy$edge)[1]),
                                                ANC = as.integer(phy$edge[,1]),
                                                DES = as.integer(phy$edge[,2])),
                         PACKAGE = "geiger")
            res=out$TIPS
            if(!length(res)) res=NULL
            return(c(node, res))
        }

    
    if(!is.rooted(tr)) stop('tree is not rooted')
    tips <- 1:length(tr$tip.label)
    intnodes <- unique(tr$edge[, 1])
    clades <- list()
    i <- 1
    while(i <= length(intnodes)){
#        descending_nodes <- get.descending.nodes.branches(tr, intnodes[i])$descending.nodes
        descending_nodes <- get.descendants.of.node(tr, intnodes[i])
        descending_tips <- tr$tip.label[descending_nodes[descending_nodes %in% tips]]
##        descending_tips <- tips(tr, intnodes[i])
        num_matches <- grep(tag, descending_tips)
        num_internal_nodes <- length(descending_nodes) - length(descending_tips)
        if(length(num_matches) == 0){
            intnodes <- intnodes[-which(intnodes %in% descending_nodes)]
        }else if(length(num_matches) == length(descending_tips)){
            clades[[length(clades)+1]] <- descending_tips
            intnodes <- intnodes[-which(intnodes %in% descending_nodes)]
        }else if(is.polytomy(tr, intnodes[i]) ){ #If it is a polytomy
            i_matches <- which(tr$tip.label %in% descending_tips[num_matches])
            ancestral_nodes <- sapply(i_matches, function(x) get.mrca(tr, x))
            tips_in_multi_clade <- i_matches[ancestral_nodes == intnodes[i]]
            if(length(tips_in_multi_clade) > 1){
                clades[[length(clades)+1]] <- tr$tip.label[tips_in_multi_clade]
            }
            intnodes <- intnodes[intnodes != intnodes[i]]
        }else{
            i <- i+1
        }
    }
    if(include.singletons){
        tips_in_clades <- unlist(clades)
        tips_out_clades <- tr$tip.label[!tr$tip.label %in% tips_in_clades]
        tips_out_clades <- grep(tag, tips_out_clades, value = T)
        return(c(clades, tips_out_clades))
    }else{
        return(clades)
    }
}

# Example
#set.seed(12344)
#tr <- rtree(10)
#tr$tip.label[c(1, 2, 3, 7, 8 , 10)]  <- paste0(tr$tip.label[c(1, 2, 3, 7, 8 , 10)], '_TAG_')

#plot(tr)
#find.monophyletic(tr, '_TAG_', include.singletons = T)
