\name{LeslieDepl}
\alias{LeslieDepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leslie Depletion Method Estimator
%%  ~~function to do ... ~~
}
\description{
Estimates abundance and probability of capture from multiple pass removal methods
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LeslieDepl(data, site.spec.p = F, Ricker.correction = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataframe or matrix with site visits across each row and counts from each pass across the columns.}
  \item{site.spec.p}{If TRUE, fits a model assuming site-specific capture probabilities as random effects. If FALSE (default), assumes common capture probability.}
  \item{Ricker.correction}{If TRUE, uses the modification proposed by Ricker in 1975, which often results in slightly higher estimates of abunance.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe consisting of 2 columns:
  \item{N.hat}{ - estimate of abundance}
  \item{N.hat.SE}{ - standard error of the estimate of abundance}
  \item{p.hat}{ - estimate of the probability of capture}

}
\references{
Ogle, D. 2011. fishR vignette - depletion methods for estimating abundance.
%% ~put references to the literature/web site here ~
}
\author{
Kevin See, QCI, Seattle, WA
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(depl.data)
# pull out examples with three-pass depletion data
three.pass.data = subset(depl.data, n.pass==3, select=c('pass.1', 'pass.2', 'pass.3'))
# estimate abundance with common catchability
LeslieDepl(three.pass.data)

# estimate abundance with site-specific catchability
LeslieDepl(three.pass.data, site.spec.p=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
