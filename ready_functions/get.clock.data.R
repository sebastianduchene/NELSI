get.clock.data <- function(rate.sim.object, tipsonly = T, ...){
  phylogram <- rate.sim.object$phylogram
  chrono <- rate.sim.object$phylogram
  chrono$edge.length <- rate.sim.object[[2]][, 7]
  times <- allnode.times(chrono, tipsonly)
  substitutions <- allnode.times(phylogram, tipsonly)
  plot(times, substitutions, ...)
  return(data.frame(times, substitutions))
}
