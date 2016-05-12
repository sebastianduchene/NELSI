dist.topo.normalised <- function(tr1, tr2, nrand=100){
    dists_sim <- vector()
    for(i in 1:nrand){
        tr_temp <- tr2
        tr_temp$tip.label <- sample(tr2$tip.label)
        dists_sim[i] <- dist.topo(tr1, tr_temp)
    }
    max_dist <- max(dists_sim)
    return(dist.topo(tr1, tr2) / max_dist)
}


