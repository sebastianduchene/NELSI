get.ltt.summary <- function(tree){
    if(is.ultrametric(tree)) stop("The tree is ultrametric. These statistics are only calculated for non-ultrametric trees")
    coords <- ltt.plot.coords(tree)
    max_location <- which.max(coords[, 2])
    relative_time_max_l <- coords[max_location, 1]/min(coords[, 1])
    reg1 <- lm(coords[1:max_location, 2] ~ coords[1:max_location, 1])
    reg2 <- lm(coords[max_location:nrow(coords), 2] ~ coords[max_location:nrow(coords), 1])
    slope1 <- reg1$coefficients[2]
    slope2 <- reg2$coefficients[2]
    ratio_slopes <- slope1/slope2
    result <- c(relative_time_max_l, slope1, slope2, ratio_slopes)
    names(result) <- c('time_max_lineages', 'slope1', 'slope2', 'ratio_slopes')
    return(result)
}