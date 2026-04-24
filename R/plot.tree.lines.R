#' Plot a phylogenetic tree using line segments
#'
#' Draws a phylogenetic tree as a series of line segments with optional
#' rotation, offset, and log-scale y axis. Unlike standard \code{plot.phylo},
#' this function uses internally computed node coordinates
#' (\code{\link{get.ordinates}}) and supports per-branch colouring.
#'
#' @param tree A rooted phylogenetic tree of class \code{"phylo"}.
#' @param rotation.angle Numeric. Angle in radians by which to rotate the
#'   entire plot. Default \code{0}.
#' @param x.offset Numeric. Horizontal shift applied after rotation. Default
#'   \code{0}.
#' @param y.offset Numeric. Vertical shift. If \code{0} (default) the minimum
#'   y coordinate is set to 0.
#' @param log.scale Logical. If \code{TRUE}, y coordinates are
#'   \code{log10(1 + y)}. Default \code{FALSE}.
#' @param line.type Character. Line type passed to \code{\link[graphics]{lines}}.
#'   Default \code{"s"} (staircase).
#' @param plot.new Logical. If \code{TRUE}, a new blank plot is opened before
#'   drawing. Default \code{FALSE}.
#' @param show.tip.labels Logical. If \code{TRUE}, tip labels are printed.
#'   Default \code{FALSE}.
#' @param lines.colour Character vector. Colours for each branch (recycled to
#'   match number of edges). Default \code{"darkgrey"}.
#' @param tips.colour Character. Colour for tip points. Default \code{"black"}.
#' @param branch.width Numeric. Line width. Default \code{2}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}} when
#'   \code{plot.new = TRUE}.
#'
#' @return \code{NULL} invisibly. Called for its plotting side effect.
#'
#' @seealso \code{\link{get.ordinates}}
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tr <- rcoal(10)
#' pdf(NULL)
#' plot.tree.lines(tr, plot.new = TRUE)
#' dev.off()
#'
#' @export plot.tree.lines
plot.tree.lines <- function(tree, rotation.angle = 0, x.offset = 0, y.offset = 0, log.scale = F,
                            line.type = "s", plot.new = F, show.tip.labels = F, lines.colour = "darkgrey",
                            tips.colour = "black", branch.width = 2, ...){
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
  if(length(lines.colour) != nrow(tree$edge)) lines.colour <- rep(lines.colour, length.out = nrow(tree$edge)) 
  for(i in 1:nrow(tree$edge)){
    ordinate_start <- ordinates_rotated$node.index == tree$edge[i, 1]
    ordinate_end <- ordinates_rotated$node.index == tree$edge[i, 2] 
    lines(c(ordinates_rotated$x.coord[ordinate_start], ordinates_rotated$x.coord[ordinate_end]),
          c(ordinates_rotated$y.coord[ordinate_start], ordinates_rotated$y.coord[ordinate_end]),
          lwd = branch.width, type = line.type, col = lines.colour[i])
  }
  tips <- ordinates_rotated$node.index %in% 1:length(tree$tip.label)
  points(ordinates_rotated$x.coord[tips], ordinates_rotated$y.coord[tips], col = tips.colour, pch = 20, cex = 0.4)
  if(show.tip.labels){
    text(ordinates_rotated$x.coord[tips], ordinates_rotated$y.coord[tips], labels = tree$tip.label, pos = 3, cex = 0.6)
  }
}
