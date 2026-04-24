#' Plot a rate simulation object
#'
#' S3 plot method for objects of class \code{"ratesim"}. Produces a two-panel
#' figure: a rate-through-time plot (one line per tip lineage) and a
#' phylogram with branch widths proportional to the log rate.
#'
#' @param rate.sim.object An object of class \code{"ratesim"} as returned by
#'   \code{\link{simulate.rate}}.
#' @param col.lineages Character vector of colours, one per tip lineage.
#'   Default uses \code{\link[grDevices]{colors}()}.
#' @param type Character. Line type passed to \code{\link[graphics]{plot}}.
#'   Default \code{"l"}.
#'
#' @return \code{NULL} invisibly. The function is called for its plotting side
#'   effect.
#'
#' @seealso \code{\link{simulate.rate}}, \code{\link{get.lineage.time.rate}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' sim <- simulate.rate(tr, simulate.clock, list(rate = 0.005, noise = 1e-6))
#' pdf(NULL)
#' plot(sim)
#' dev.off()
#'
#' @export
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
