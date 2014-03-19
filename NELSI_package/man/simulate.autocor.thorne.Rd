\name{simulate.autocor.thorne}
\alias{simulate.autocor.thorne}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulate.autocor.thorne(tree, params = list(initial.rate = 0.01, v = 0.3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tree, params = list(initial.rate = 0.01, v = 0.3)) 
{
    require(phangorn)
    require(geiger)
    initial.rate <- params$initial.rate
    v = params$v
    data.matrix <- get.tree.data.matrix(tree)
    while (any(is.na(data.matrix[, 5])) | any(is.nan(data.matrix[, 
        5]))) {
        data.matrix[1, 5] <- initial.rate
        for (i in 2:nrow(data.matrix)) {
            parent.node <- data.matrix[i, 2]
            preceeding.parent <- data.matrix[, 2][data.matrix[, 
                3] == parent.node]
            preceeding.parent.brage <- data.matrix[, 4][data.matrix[, 
                2] == preceeding.parent][1]
            preceeding.parent.brrate <- data.matrix[, 5][data.matrix[, 
                2] == preceeding.parent][1]
            if (!(is.na(preceeding.parent.brrate)) & !(is.nan(preceeding.parent.brrate)) & 
                (parent.node \%in\% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(preceeding.parent.brrate)), 
                  sd = v * abs(data.matrix[i, 4] - preceeding.parent.brage)^0.5))
            }
            else if (!(parent.node \%in\% data.matrix[, 3])) {
                data.matrix[i, 5] <- abs(rlnorm(1, mean = log(abs(initial.rate)), 
                  sd = sqrt(initial.rate)))
            }
        }
    }
    data.matrix[, 6] <- data.matrix[, 7] * data.matrix[, 5]
    tree$edge.length <- data.matrix[, 6]
    res <- list(tree, data.matrix)
    names(res) <- c("phylogram", "tree.data.matrix")
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
