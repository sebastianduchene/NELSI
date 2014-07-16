
library(NELSI)
library(phangorn)
source('incomplete_trees.R')

# Define all simulation settings:

sc_rate <- 0.005
sc_noise <- 0.0001

ucl_mean <- log(0.005)
ucl_sd <- 0.1

uch_mean <- log(0.005)
uch_sd <- 0.3

ucgl_alpha <- 10
ucgl_rate <- 2000

ucgh_alpha <- 1
ucgh_rate <- 200

acl_init <- 0.005
acl_v <- 0.3

acm_init <- 0.005
acm_v <- 0.2

ach_init <- 0.005
ach_v <- 0.1


trees_names <- paste0("complete_trees/comp_yule_", 1:10, ".tre")

trees_phy <- list()

# Read trees
for(i in 1:length(trees_names)){
  trees_phy[[i]] <- read.tree(trees_names[i])
}
names(trees_phy) <- gsub('complete_trees/', '',trees_names)

# Find taxa to remove from each tree
taxa_remove <- list()
for(i in 1:length(trees_phy)){
  taxa_remove[[i]] <- prune_tree_random(trees_phy[[i]], 50)
  cat(paste0(taxa_remove[[i]], collapse = ' '), file = 'taxa_remove.txt', sep = '\n', append = T)
}
names(taxa_remove) <- names(trees_phy)



for(i in 1:length(trees_phy)){
  
  # Start strict clock
  print(paste("simulating strict clock on", names(trees_phy)[i]))
  var_sites <- 0
  counter_i <- 0
  counter_tot <- 0   
  while(var_sites > 5200 || var_sites < 4800){
  break
    print(paste("rep", counter_i, "params are", sc_rate, var_sites))
    sim_tree <- simulate.clock(trees_phy[[i]] , params = list(rate = sc_rate, noise = sc_noise))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 10){
      if(var_sites > 5200){
        sc_rate <- sc_rate - (sc_rate*0.2)
        counter_i <- 0
      }else if(var_sites < 4800){
          sc_rate <- sc_rate + (sc_rate*0.2)
	  counter_i <- 0
      }
    }
  }

  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('sc_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  # End strict clock

  # Start UCL
    print(paste("simulating UCL clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){
  break
    print(paste("rep", counter_i, "params are", ucl_mean, var_sites))
    sim_tree <- simulate.uncor.lnorm(trees_phy[[i]] , params = list(mean.log = ucl_mean, sd.log = ucl_sd))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        ucl_mean <- log(exp(ucl_mean) - (exp(ucl_mean)*0.2))
        counter_i <- 0
      }else if(var_sites < 4800){
          ucl_mean <- log(exp(ucl_mean) + (exp(ucl_mean)*0.2))
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('ucl_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End UCL

  # Start UCH
    print(paste("simulating UCH clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){
break  
    print(paste("rep", counter_i, "params are", uch_mean, var_sites))
    sim_tree <- simulate.uncor.lnorm(trees_phy[[i]] , params = list(mean.log = uch_mean, sd.log = uch_sd))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        uch_mean <- log(exp(uch_mean) - (exp(uch_mean)*0.2))
        counter_i <- 0
      }else if(var_sites < 4800){
          uch_mean <- log(exp(uch_mean) + (exp(uch_mean)*0.2))
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('uch_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End UCH


  # Start UCGL
    print(paste("simulating UCGL clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){
break
    sim_tree <- simulate.uncor.gamma(trees_phy[[i]] , params = list(shape = ucgl_alpha, rate = ucgl_rate))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    print(paste("rep", counter_i, "params are", ucgl_rate, var_sites))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        ucgl_rate <- ucgl_rate + (ucgl_rate*0.1)
        counter_i <- 0
      }else if(var_sites < 4800){
          ucgl_rate <- ucgl_rate - (ucgl_rate*0.1)
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('ucgl_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End UCGL



  # Start UCGH
    print(paste("simulating UCGH clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){
break
    sim_tree <- simulate.uncor.gamma(trees_phy[[i]] , params = list(shape = ucgh_alpha, rate = ucgh_rate))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    print(paste("rep", counter_i, "params are", ucgh_rate, var_sites))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        ucgh_rate <- ucgh_rate + (ucgh_rate*0.1)
        counter_i <- 0
      }else if(var_sites < 4800){
          ucgh_rate <- ucgh_rate - (ucgh_rate*0.1)
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('ucgh_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End UCGH



  # Start ACL
    print(paste("simulating ACL clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){
break
    sim_tree <- simulate.autocor.kishino(trees_phy[[i]] , params = list(initial.rate = acl_init, v = acl_v))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    print(paste("rep", counter_i, "params are", acl_init, var_sites))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        acl_init <- acl_init - (acl_init*0.1)
        counter_i <- 0
      }else if(var_sites < 4800){
          acl_init <- acl_init + (acl_init*0.1)
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('acl_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End ACL


  # Start ACM
    print(paste("simulating ACM clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){
break
    sim_tree <- simulate.autocor.kishino(trees_phy[[i]] , params = list(initial.rate = acm_init, v = acm_v))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    print(paste("rep", counter_i, "params are", acm_init, var_sites))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        acm_init <- acm_init - (acm_init*0.1)
        counter_i <- 0
      }else if(var_sites < 4800){
          acm_init <- acm_init + (acm_init*0.1)
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('acm_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End ACM


# TESTED UP TO HERE
  # Start ACH
    print(paste("simulating ACH clock on", names(trees_phy)[i]))
    var_sites <- 0
    counter_i <- 0
    counter_tot <- 0
  while(var_sites > 5200 || var_sites < 4800){

    sim_tree <- simulate.autocor.kishino(trees_phy[[i]] , params = list(initial.rate = ach_init, v = ach_v))$phylogram
    sim_prune <- ape::drop.tip(sim_tree, taxa_remove[[i]])
    sim_data <- as.DNAbin(simSeq(sim_prune, l = 10000))
    var_sites <- length(seg.sites(sim_data))
    print(paste("rep", counter_i, "params are", ach_init, var_sites))
    counter_tot <- counter_tot + 1
    if(counter_tot > 200){
      print("aborting simulation")
      break
    }
    counter_i <- counter_i + 1
    if(counter_i > 5){
      if(var_sites > 5200){
        ach_init <- ach_init - (ach_init*0.1)
        counter_i <- 0
      }else if(var_sites < 4800){
          ach_init <- ach_init + (ach_init*0.1)
          counter_i <- 0
      }
    }
  }
    
  if(var_sites < 5200 && var_sites > 4800){
    write.dna(sim_data, file= paste0('ach_', gsub('tre', 'fasta', names(trees_phy)[i])), format = "fasta", nbcol = -1, colsep = '')
  }
  #End ACH



}
