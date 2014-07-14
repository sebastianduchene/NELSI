library(phangorn)
sim_trees <- grep('yule*', dir(), value = T)

for(i in sim_trees){
      tr_temp <- read.tree(i)
      tip_names <- tr_temp$tip.label
      if("outgroup" %in% tip_names){
        print(paste("tree", i, "has an outgroup. I will remove it"))
        tr_temp <- drop.tip(phy = tr_temp, tip = 'outgroup')
	write.tree(tr_temp, file = i)
      }else{
        print(paste("tree", i, "has no outgroup. Moving to next"))
	next
      }
}
      