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
  \item{rmInvalid}{if TRUE, marks any invalid estimates of N.hat and N.hat.SE as NA}
}
\details{
The data must consiste of a dataframe or matrix with three columns. The first is the number of individuals marked after the first capture event. The second is the total number of individuals captured during the second capture (recapture) event. The third is the number of those captured during the second capture event that were marked individuals, i.e. the number of recaptures.

The equation used is:
N = [(M+1) * (C+1)] / (R+1) - 1

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
data(MR.data)
# estimate abundance
ChapmanMR(MR.data[,c('M', 'C', 'R')])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
