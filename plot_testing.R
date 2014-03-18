setwd("ready_functions")
for(i in dir()) source(i)
setwd("..")

library(phangorn)
library(geiger)

tr1 <- sim.bdtree(b = 1, d = 0, n = 10, stop = "taxa")

conv.factor <- max(branching.times(tr1)) / 50

tr1$edge.length <- tr1$edge.length / conv.factor

r1 <- simulate.autocor.kishino(tr1)



rates.time.list <- list()
for(i in 1:length(r1[[1]])){
      rates.time.list[[i]]<- get.lineage.time.rate(i, r1[[2]])
}

ylims <- range(sapply(rates.time.list, function(y) range(y[,2])))

chrono <- r1[[1]]
chrono$edge.length <- r1[[2]][, 7]

xlims <- range(branching.times(chrono$edge.length))





#plot(rates.time.list[[1]], type = "l", col = rgb(0, 0, 1, 0.5), lwd = 3)


#plot.rates.time <- function(


