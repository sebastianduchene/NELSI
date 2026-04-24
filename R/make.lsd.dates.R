#' Write a sampling-date file for LSD
#'
#' Extracts sampling dates from tip labels and writes them to a tab-delimited
#' file in the format expected by LSD (Least Squares Dating). Dates are
#' extracted by stripping the portion of each tip label matched by
#' \code{grep.sep}.
#'
#' @param tr A phylogenetic tree of class \code{"phylo"} whose tip labels
#'   encode sampling dates (e.g. \code{"taxon_2015.5"}).
#' @param grep.sep Character. Regular expression matched and removed from each
#'   tip label to obtain the date string. Default \code{".+_"} strips
#'   everything up to and including the last underscore.
#' @param outfile Character. Path to the output file. Default \code{"out.date"}.
#' @param random Logical. If \code{TRUE}, dates are randomly shuffled across
#'   tips before writing. Default \code{FALSE}.
#'
#' @return \code{NULL} invisibly. The function is called for its side effect
#'   of writing \code{outfile}.
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(5)
#' tr$tip.label <- paste0("taxon_", c(2010, 2011, 2012, 2013, 2014))
#' make.lsd.dates(tr, outfile = tempfile())
#'
#' @export
make.lsd.dates <- function(tr, grep.sep = '.+_', outfile = 'out.date', random = F){
  dates <- gsub(grep.sep, '', tr$tip.label)
  if(random){
    dates <- sample(dates)
  }
  dates.vector <- paste(tr$tip.label, dates, sep = '\t')
  cat(length(tr$tip.label), '\n', file = outfile)
  cat(dates.vector, sep = '\n', file = outfile, append = T)
}
