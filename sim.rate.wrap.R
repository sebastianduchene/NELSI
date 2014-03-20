# Tesing a function that wrapps all the rate simulation functions. It should return an object of class ratesim

setwd("ready_functions")
for(i in dir()) source(i)
setwd("..")

library(phangorn)


tree <- rcoal(10)
params <- list(initial.rate = 0.01, v = 0.3)


r1 <- simulate.autocor.thorne(tree = tree, params = params)


#Function starts here. 

simulate.rate <- function(tree, FUN, ...){
	 ratesim.object <- FUN(tree, ...)
	 return(ratesim.object)
}



#expe1

e1 <- sim.rate(tree, FUN = simulate.uncor.lnorm, params = list(mean.log = -3.9, sd.log = 0.8))

plot(e1)

system("sleep 2")

#expe 2

e2 <- sim.rate(tree, FUN = simulate.clock, params = list(rate = 0.02, noise = 0))

plot(e2)