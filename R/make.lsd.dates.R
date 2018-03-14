make.lsd.dates <- function(tr, grep.sep = '.+_', outfile = 'out.date', random = F){
  dates <- gsub(grep.sep, '', tr$tip.label)
  if(random){
    dates <- sample(dates)
  }
  dates.vector <- paste(tr$tip.label, dates, sep = '\t')
  cat(length(tr$tip.label), '\n', file = outfile)
  cat(dates.vector, sep = '\n', file = outfile, append = T)
}
