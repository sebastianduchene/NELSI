

library(phangorn)

data_fil <- grep('fasta', dir(), value = T)

for(i in 1:length(data_fil)){

print(paste('testing', data_fil[i]))

data_1 <- read.dna(data_fil[i], format = "fasta")

n_tree <- as.numeric(gsub('[A-Z]|[a-z]| |[.]|_', '', data_fil[i]))

print(n_tree)


tree_comp_1 <- read.tree(paste0('../complete_trees/comp_yule_', n_tree, '.tre'))


tree_prune_1 <- read.tree(paste0('../incomplete_trees/pruned_yule_', n_tree, '.tre'))


# Which are the pruned taxa?

pruned_taxa <- tree_comp_1$tip.label[!(tree_comp_1$tip.label %in% tree_prune_1$tip.label)]

# Which taxa are kept?

q1 <- all(tree_comp_1[tree_comp_1$tip.label %in% tree_prune_1] %in% rownames(data_1) )
print(q1)

kept_taxa <- tree_prune_1$tip.label

# Are all the kept taxa in the dna matirx

q2 <- all(rownames(data_1) %in% kept_taxa)

print(q2)

if(!(q1 && q2)) stop(paste('data set', data_fil[i], 'has a taxon mismatch problem'))

#print(pruned_taxa)

}