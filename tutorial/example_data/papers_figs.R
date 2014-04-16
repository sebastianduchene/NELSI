#uncor <- simulate.uncor.lnorm(tr, params = list(mean.log = 0.03, sd.log = 1))

#kishino <- simulate.autocor.kishino(tr, params = list(initial.rate = 0.02, v = 0.5))

#clock <- simulate.clock(tree = tr, params = list(rate = 0.02, noise = 0))

par(mfrow = c(1, 3))
plot(tr, show.tip.label = F, edge.width = log10(clock$tree.data.matrix[, 5]/0.00001))

plot(tr, show.tip.label = F, edge.width = log10(kishino$tree.data.matrix[, 5]/0.000001 ))

plot(tr, show.tip.label = F, edge.width = log10(uncor$tree.data.matrix[, 5]/0.0001 ))


dev.copy2pdf()


plot(clock$phylogram, show.tip.label = F, edge.width = 2)
plot(kishino$phylogram, show.tip.label = F, edge.width = 2)
plot(uncor$phylogram, show.tip.label = F, edge.width = 2)
