plot.ratesim <-
function(rate.sim.object, col.lineages = colors(), type = "l"){
  
  rates.time.list <- list()
  for(i in 1:length(rate.sim.object[[1]]$tip.label)){
      rates.time.list[[i]]<- get.lineage.time.rate(i, rate.sim.object)
   }

   ylims <- range(lapply(rates.time.list, function(y) range(y[,2])))
   chrono <- rate.sim.object[[1]]
   chrono$edge.length <- rate.sim.object[[2]][, 7]
   node.ages <- allnode.times(chrono)
   xlims <- sort(range(node.ages), decreasing = T)      

   par(mfrow = c(1, 2))
   plot(rates.time.list[[1]][, 1], rates.time.list[[1]][, 2], ylim = ylims, xlim = xlims,  ylab = "Rate", xlab = "Time", type = type, lwd = 3, col = col.lineages[1])
   for(k in 2:length(rates.time.list)){
     lines(rates.time.list[[k]][, 1], rates.time.list[[k]][, 2], ylim = ylims, xlim = xlims, col = col.lineages[k], lwd = 3, type = type) 
     }

     plot(chrono, edge.width = 1 + log( rate.sim.object[[2]][, 5] / min(rate.sim.object[[2]][, 5])), show.tip.label = F, root.edge = T)
     tiplabels(pch = 16, col = col.lineages, cex = 1.5)
}
