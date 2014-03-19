
library(phangorn)
library(geiger)

t1 <- sim.bdtree(b = 1, d = 0, n = 10, stop = "taxa")

setwd("ready_functions")
for(i in dir()) source(i)
setwd("..")

r1 <- simulate.tdep.generic(t1, params = list(mu = 0.015, srate = 0.035, lambda = 0.2, noise = 0.0002))

plot.rates.time(rate.sim.object = r1, col.lineages = rainbow(10))
