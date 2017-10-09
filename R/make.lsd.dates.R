make.lsd.dates <- function(tr, youngest_sample = 0, output_file){
  tip_ages <- abs(allnode.times(tr, tipsonly = T) - youngest_sample)
  ages_lines <- c(length(tip_ages), paste(tr$tip.label, tip_ages, sep = '\t'))
  writeLines(ages_lines, con = output_file)
}