

setwd("ready_functions")
for(i in dir()) source(i)
setwd("..")

library(geiger)

tr1 <- rtree(10)

rate.sim.object <- simulate.autocor.kishino(tr1, params = list(initial.rate = 0.01, v = 0.5))


plot.clock <- function(rate.sim.object, tipsonly = T, ...){
  phylogram <- rate.sim.object$phylogram
  chrono <- rate.sim.object$phylogram
  chrono$edge.length <- rate.sim.object[[2]][, 7]
  times <- allnode.times(chrono, tipsonly)
  substitutions <- allnode.times(phylogram, tipsonly)
  plot(times, substitutions, ...)
  return(data.frame(times, substitutions))
}


test1 <- plot.clock(rate.sim.object, pch = 20, col = "red")
