\name{make_lsd_date_file}
\alias{make_lsd_date_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make_lsd_date_file
}
\description{
make_lsd_date_file makes a .date file for LSD

}
\usage{
make_lsd_date_file(phylodata, outfile)
}
\arguments{
  \item{phylodata}{
A phylogenetic tree of class phylo or a sequence alignment of class DNAbin
}
}
\details{
Does not return anything. Instead, it saves the dates in the files specified through outfile
}
\value{
None.
}
\references{
Pending.
}
\author{
Sebastian Duchene
}
\note{
None
}

\seealso{

}
\examples{
Pending.

## The function is currently defined as

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


}
\keyword{ chronogram }
