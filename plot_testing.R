setwd("ready_functions")
for(i in dir()) source(i)
setwd("..")

library(phangorn)
library(geiger)

tr1 <- sim.bdtree(b = 1, d = 0, n = 10, stop = "taxa")

conv.factor <- max(branching.times(tr1)) / 50

tr1$edge.length <- tr1$edge.length / conv.factor

r1 <- simulate.autocor.kishino(tr1, params = list(initial.rate = 0.01, v = 0.1))


col.lineages <- rgb(0, 0, 1, 0.5)

col.lineages <- rainbow(length(r1[[1]]$tip.label))

col.lineages <- terrain.colors(length(r1[[1]]$tip.label))

col.lineages <- heat.colors(length(r1[[1]]$tip.label))

col.lineages <- sample(size = length(r1[[1]]$tip.label), colors()) 

type = "l"

# Function starts here



rates.time.list <- list()
for(i in 1:length(r1[[1]]$tip.label)){
      rates.time.list[[i]]<- get.lineage.time.rate(i, r1[[2]])
}

ylims <- range(sapply(rates.time.list, function(y) range(y[,2])))

chrono <- r1[[1]]
chrono$edge.length <- r1[[2]][, 7]

xlims <- sort(range(branching.times(chrono)), decreasing = T)

par(mfrow = c(1, 2))
plot(rates.time.list[[1]][, 1], rates.time.list[[1]][, 2], ylim = ylims, xlim = xlims,  ylab = "Rate", xlab = "Time", type = type, lwd = 3)

for(k in 2:length(rates.time.list)){
     lines(rates.time.list[[k]][, 1], rates.time.list[[k]][, 2], ylim = ylims, xlim = xlims, col = col.lineages[k], lwd = 3, type = type) 
}

plot(chrono, edge.width = 1 + log( r1[[2]][, 5] / min(r1[[2]][, 5])), show.tip.label = F)#, tip.color = col.lineages)
tiplabels(pch = 16, col = col.lineages, cex = 1.5)


#plot(rates.time.list[[1]], type = "l", col = rgb(0, 0, 1, 0.5), lwd = 3)


#plot.rates.time <- function(


