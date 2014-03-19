\name{get.rate.descendant.pairs}
\alias{get.rate.descendant.pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.rate.descendant.pairs(tree.data.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree.data.matrix}{
%%     ~~Describe \code{tree.data.matrix} here~~
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
function (tree.data.matrix) 
{
    dat <- tree.data.matrix
    parent.rate <- vector()
    daughter.rate <- vector()
    br.len <- vector()
    for (i in 1:nrow(dat)) {
        daughter.temp <- dat[i, 3]
        parent.temp <- dat[i, 2]
        br.temp <- dat[i, 4]
        if (parent.temp \%in\% dat[, 3]) {
            parent.rate <- c(parent.rate, dat[dat[, 3] == parent.temp, 
                5])
            daughter.rate <- c(daughter.rate, dat[i, 5])
            br.len <- c(br.len, abs(br.temp - dat[dat[, 3] == 
                parent.temp, 5]))
        }
    }
    return(data.frame(parent.rate, daughter.rate, br.len))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
