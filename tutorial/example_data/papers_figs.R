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

hivTree <- read.annotated.nexus("hiv_A_env.tree")

par(mar = c(4, 3.8, 4, 3.5))
par(mfrow = c(1, 2))
plot(hivTree, cex = 0.7, label.offset = 7)
axisPhylo()


nodeAges <- allnode.times(hivTree)
tiplabels(round(nodeAges[1:length(hivTree$tip.label)], 2), cex = 0.7, adj = c(0.25, 0.5))
nodelabels(round(nodeAges[(length(hivTree$tip.label) +1):length(nodeAges)], 0), cex = 0.8)

hivData <- trann2trdat(hivTree)

hivPhylogram <- hivTree # copy the the HIV tree into an other variable
hivPhylogram$edge.length <- hivData$blensubs # use the substitutions for the branch lengths
nodeTipSubst <- allnode.times(hivPhylogram, tipsonly = T)
nodeTipTime <- allnode.times(hivTree, tipsonly = T)


plot(nodeTipTime, nodeTipSubst, pch = 20, ylab = "Substitutions along a lineage (subst)", xlab = "Time along a lineage (years)")
hivLM <- lm(nodeTipSubst ~ nodeTipTime)
abline(hivLM)
summary(hivLM)

