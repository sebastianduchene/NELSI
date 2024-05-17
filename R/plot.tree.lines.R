plot.tree.lines <- function(tree, rotation.angle = 0, x.offset = 0, y.offset = 0, log.scale = F, 
                            line.type = "s", plot.new = F, show.tip.labels = F, lines.colour = "darkgrey",
                            points.colour = "black", ...){
    rotate <- function(v, angle){
        tmatrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2, byrow=T)
        return(tmatrix %*% v)
    }
  ordinates <- get.ordinates(tree)
  ordinates_rotated <- ordinates
  for(i in 1:nrow(ordinates_rotated)){
    ordinates_rotated[i, c("x.coord", "y.coord")] <- rotate(c(ordinates_rotated[i, "x.coord"], ordinates_rotated[i, "y.coord"]), rotation.angle)
  }
  ordinates_rotated <- data.frame(ordinates_rotated)
  if(y.offset == 0){
    ordinates_rotated$y.coord <- ordinates_rotated$y.coord - min(ordinates_rotated$y.coord)
  } else {
    ordinates_rotated$y.coord <- ordinates_rotated$y.coord + y.offset
  }
  if(x.offset != 0){
    ordinates_rotated$x.coord <- ordinates_rotated$x.coord + x.offset
  }
  if(log.scale){
    ordinates_rotated$y.coord <- log10(1 + ordinates_rotated$y.coord)
  }
  if(plot.new){
    plot(c(min(ordinates_rotated$x.coord), max(ordinates_rotated$x.coord)), 
         c(min(ordinates_rotated$y.coord), max(ordinates_rotated$y.coord)), 
         type = "n", bty = "n", ...)
  }
  for(i in 1:nrow(tree$edge)){
    ordinate_start <- ordinates_rotated$node.index == tree$edge[i, 1]
    ordinate_end <- ordinates_rotated$node.index == tree$edge[i, 2]
    lines(c(ordinates_rotated$x.coord[ordinate_start], ordinates_rotated$x.coord[ordinate_end]),
          c(ordinates_rotated$y.coord[ordinate_start], ordinates_rotated$y.coord[ordinate_end]),
          col = "darkgrey", lwd = 2, type = line.type, col = lines.colour)
  }
  tips <- ordinates_rotated$node.index %in% 1:length(tree$tip.label)
  points(ordinates_rotated$x.coord[tips], ordinates_rotated$y.coord[tips], col = tips.colour, pch = 20, cex = 0.4)
  if(show.tip.labels){
    text(ordinates_rotated$x.coord[tips], ordinates_rotated$y.coord[tips], labels = tree$tip.label, pos = 3, cex = 0.6)
  }
}
