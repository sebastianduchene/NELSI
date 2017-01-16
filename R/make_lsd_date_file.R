make_lsd_date_file <- function(phylodata, outfile = 'outfile.date'){
    if(class(phylodata) == 'DNAbin'){
        taxa_names <- rownames(phylodata)
    }else if(class(phylodata) == 'phylo'){
        taxa_names <- phylodata$tip.label
    }

    dates <- sapply(taxa_names, function(x) gsub('.+_', '', x), USE.NAMES = F)
    lines <- paste0(taxa_names, ' ', dates, collapse = '\n')
    cat(length(taxa_names), '\n', file = outfile)
    cat(lines, file = outfile, append = T)
    print(paste('Dates file saved in ', outfile))
}
