
print("
The original simulations are:

    - SC
      rate = 0.005, noise = rnorm(n, 0, 0.0001)
    - UCLN-L
      mean.rate = 0.005, sd = 0.1
    - UCLN-H
      mean.rate = 0.005, sd = 0.3
    - UCG-L # shape = alpha, rate = beta
      alpha = 10, beta = 10^5 (10^4)
    - UCG-H
      alpha = 1, beta = 10^4 (10^3)
    - ACLN-L
      rate = 0.005, v = 0.3
    - ACLN-M
      rate = 0.005, v = 0.2
    - ACLN-H
      rate = 0.005, v = 0.1

Ten data sets per simulation settings, with a root age of 10 time units. In all cases, the number of variable sites is between 5300 and 4800 in alignments of 10000 variable sites

Use rates that are 1/5 of magnitude lower, as Reviewer 1 suggests.
#    TODO   
    - Check sequence variation. It should be around 1000 variables sites
    - Simulate yule trees with 100 tips, and 10 time units. 
")

#eval(parse(text = "print(tr1)")


#################
#################
library(phangorn)

cudir <- getwd()
setwd("~/Desktop/nelsi/NELSI/R")
for(i in dir()) source(i)
setwd(cudir)

print("Running ucl simulations")

ucl_rate <- 0.001
ucl_sd <- 0.1


tree_name <- 'yule1.tre'
tree_temp <- read.tree(tree_name)
tree_temp$edge.length <- tree_temp$edge.length * (10 / max(branching.times(tree_temp)))

phylogram_temp <- tree_temp
sim_ucl <- simulate.uncor.lnorm(tree_temp, params = list(mean.log = ucl_rate, sd.log = 0.1))
## pending



stop('wow')

 

var_sites <- 0
i <- 1
while(var_sites > 1200 || var_sites < 850 ){
  phylogram_temp$edge.length <- abs((tree_temp$edge.length * sc_rate) + rnorm(length(phylogram_temp$edge.length), 0, 0.00001))
  seq_data <- as.DNAbin(simSeq(phylogram_temp, l = 10000))
  var_sites <- length(seg.sites(seq_data))
  print(paste("The number of variable sites is", var_sites))
  print(paste("Generating data with specified variables sites. Replicate ", i))
  i <- i + 1
  if(i > 50){
    print("reached maximum replicates")
    break
  }
}

if(var_sites < 1200 || var_sites > 850){
  write.dna(seq_data, file = gsub('tre', 'fasta', paste0('sc_', tree_name)), format = 'fasta', nbcol = -1, colsep = '')
}