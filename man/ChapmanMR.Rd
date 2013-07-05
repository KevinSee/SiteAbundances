\name{ChapmanMR}
\alias{ChapmanMR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chapman Mark/Recapture Abundance Estimator
%%  ~~function to do ... ~~
}
\description{
This is the modified Lincoln-Petersen or Chapman estimator. It estimates abundance based on a single mark event with a single recapture event.
}
\usage{
ChapmanMR(data, rmInvalid = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe or matrix with 3 columns for the number of marked, captured and recaptured individuals.}
  \item{rmInvalid}{if TRUE, removes any invalid estimates of N.hat and N.hat.SE}
}
\details{
The data must consiste of a dataframe or matrix with three columns. The first is the number of individuals marked after the first capture event. The second is the total number of individuals captured during the second capture (recapture) event. The third is the number of those captured during the second capture event that were marked individuals, i.e. the number of recaptures.

Estimates may be considered invalid when the product of the number of individuals captured on the first and second capture events is smaller than 4 times the estimated abundance. Based on Robson & Regier criteria.
}
\value{
    \item{N.hat}{estimate of abundance}
    \item{N.hat.SE}{standard error of abundance estimate}
    \item{p.hat}{estimate of the probability of capture}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kevin See, QCI, Seattle, WA
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
function (data, rmInvalid = T) 
{
    if (is.data.frame(data) == F) 
        data = as.data.frame(data)
    names(data) = c("M", "C", "R")
    N.hat = with(data, ((M + 1) * (C + 1))/(R + 1) - 1)
    N.hat.SE = with(data, sqrt(((M + 1) * (C + 1) * (M - R) * 
        (C - R))/((R + 1)^2 * (R + 2))))
    p.hat = with(data, R/M)
    Valid = with(data, (M * C)) > (N.hat * 4)
    if (rmInvalid == T) {
        N.hat[Valid == F] = NA
        N.hat.SE[Valid == F] = NA
    }
    return(data.frame(N.hat = N.hat, N.hat.SE = N.hat.SE, p.hat = p.hat, 
        Valid = Valid))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
