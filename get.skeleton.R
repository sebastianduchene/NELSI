library(geiger)
library(phangorn)

setwd("ready_functions")
for(i in dir()) source(i)
setwd("..")

package.skeleton(name = "NELSI_0.1")